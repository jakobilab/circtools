#!/usr/bin/env python3

# Copyright (C) 2017 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Optional, anonymous first-run usage survey.

On the first interactive run, the user is asked whether they want to share
four pieces of information: operating system, OS version, circtools version
and Python version. Nothing else is collected. The answer (yes or no) is
stored locally so the question is only ever asked once.

Design rules:
  * Uses only the standard library (no new dependencies for circtools).
  * Never blocks or breaks the real command: every failure is swallowed.
  * Never prompts in non-interactive contexts (pipes, cron, CI, containers
    without a TTY, HPC batch jobs).
  * Opt out entirely with the environment variable CIRCTOOLS_NO_SURVEY=1.
  * Override the endpoint with CIRCTOOLS_SURVEY_URL (useful for testing).
"""

import json
import os
import platform
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# TODO: replace with the real production endpoint once the server is deployed
DEFAULT_SURVEY_URL = "http://10.228.21.77:8111/api/v1/survey"

OPT_OUT_ENV = "CIRCTOOLS_NO_SURVEY"
URL_ENV = "CIRCTOOLS_SURVEY_URL"

SCHEMA_VERSION = 1
REQUEST_TIMEOUT_SECONDS = 3

# Arguments for which we never want to interrupt the user
_PASSTHROUGH_ARGS = {"-h", "--help", "-V", "--version"}


# --------------------------------------------------------------------------- #
# local state
# --------------------------------------------------------------------------- #

def _state_file():
    config_home = os.environ.get("XDG_CONFIG_HOME")
    base = Path(config_home) if config_home else Path.home() / ".config"
    return base / "circtools" / "survey.json"


def _load_state():
    try:
        with open(_state_file(), "r") as handle:
            state = json.load(handle)
        return state if isinstance(state, dict) else None
    except (OSError, ValueError):
        return None


def _save_state(state):
    """Write state atomically. Returns True on success."""
    path = _state_file()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        with open(tmp, "w") as handle:
            json.dump(state, handle, indent=2)
        os.replace(tmp, path)
        return True
    except OSError:
        return False


# --------------------------------------------------------------------------- #
# data collection
# --------------------------------------------------------------------------- #

def _os_info():
    """Return (os_name, os_version) in a human-friendly form."""
    system = platform.system()

    if system == "Darwin":
        return "macOS", platform.mac_ver()[0] or platform.release()

    if system == "Linux":
        # Prefer the distribution name (e.g. Ubuntu 22.04) over the kernel
        try:
            release = platform.freedesktop_os_release()  # Python 3.10+
            name = release.get("NAME") or "Linux"
            ver = release.get("VERSION_ID") or release.get("BUILD_ID") or ""
            if ver:
                return name, ver
        except (AttributeError, OSError):
            pass
        return "Linux", platform.release()

    if system == "Windows":
        return "Windows", platform.version()

    return system or "unknown", platform.release() or "unknown"


def collect_survey_data(circtools_version):
    os_name, os_version = _os_info()
    return {
        "schema_version": SCHEMA_VERSION,
        "os": os_name,
        "os_version": os_version,
        "circtools_version": str(circtools_version),
        "python_version": platform.python_version(),
    }


# --------------------------------------------------------------------------- #
# network
# --------------------------------------------------------------------------- #

def _send(data):
    """POST the survey data. Returns True only on a 2xx response."""
    url = os.environ.get(URL_ENV, DEFAULT_SURVEY_URL)
    request = urllib.request.Request(
        url,
        data=json.dumps(data).encode("utf-8"),
        headers={"Content-Type": "application/json",
                 "User-Agent": "circtools-survey/%d" % SCHEMA_VERSION},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            return 200 <= response.status < 300
    except (urllib.error.URLError, OSError, ValueError):
        return False


# --------------------------------------------------------------------------- #
# user interaction
# --------------------------------------------------------------------------- #

def _is_interactive():
    try:
        return sys.stdin.isatty() and sys.stderr.isatty()
    except (AttributeError, ValueError):
        return False


def _ask_user(data):
    """Show what would be sent and ask for consent.

    Returns True (yes), False (no) or None (no answer, e.g. Ctrl-C/EOF).
    """
    err = sys.stderr
    err.write(
        "\n"
        "circtools would like to ask for a small favor (this is shown only once).\n"
        "\n"
        "Would you like to take part in a short, anonymous usage survey? It sends\n"
        "only the following information to the circtools developers, so we know\n"
        "which platforms to support:\n\n"
        "    Operating system:    %s\n"
        "    OS version:          %s\n"
        "    circtools version:   %s\n"
        "    Python version:      %s\n"
        "\n"
        "Nothing else is collected. To disable this question permanently without\n"
        "answering, set the environment variable %s=1.\n"
        "\n" % (data["os"], data["os_version"], data["circtools_version"],
                data["python_version"], OPT_OUT_ENV)
    )
    err.flush()

    while True:
        try:
            err.write("Participate in the survey? [y/N]: ")
            err.flush()
            answer = input().strip().lower()
        except (EOFError, KeyboardInterrupt):
            err.write("\n")
            return None

        if answer in ("y", "yes"):
            return True
        if answer in ("", "n", "no"):
            return False
        err.write("Please answer 'y' or 'n'.\n")


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #

def maybe_run_survey(circtools_version, argv=None):
    """Call once at start-up. Safe to call on every run; never raises."""
    try:
        _maybe_run_survey(circtools_version, argv if argv is not None else sys.argv[1:])
    except Exception:
        # The survey must never get in the way of the actual analysis
        pass


def _maybe_run_survey(circtools_version, argv):
    if os.environ.get(OPT_OUT_ENV):
        return

    # don't interrupt `circtools --help`, `circtools -V` or a bare `circtools`
    if not argv or argv[0] in _PASSTHROUGH_ARGS:
        return

    state = _load_state()

    if state is not None:
        # already answered: declined, or accepted and delivered -> nothing to do
        if state.get("status") == "declined" or state.get("sent"):
            return

        # accepted earlier but the upload failed (e.g. offline): retry quietly
        if state.get("status") == "accepted":
            if _send(collect_survey_data(circtools_version)):
                state["sent"] = True
                state["sent_at"] = int(time.time())
                _save_state(state)
            return

    # first run: only ask if a human can actually answer
    if not _is_interactive():
        return

    # if we can't persist the answer we would nag on every run, so don't ask
    if not _save_state({"status": "pending", "asked_at": int(time.time())}):
        return

    data = collect_survey_data(circtools_version)
    answer = _ask_user(data)

    if answer is None:
        # Ctrl-C / EOF: not an answer, ask again next time
        try:
            _state_file().unlink()
        except OSError:
            pass
        return

    if answer is False:
        _save_state({"status": "declined", "answered_at": int(time.time())})
        sys.stderr.write("No problem, you won't be asked again.\n\n")
        return

    state = {"status": "accepted", "answered_at": int(time.time()), "sent": False}
    if _send(data):
        state["sent"] = True
        state["sent_at"] = int(time.time())
        sys.stderr.write("Thank you for participating!\n\n")
    else:
        sys.stderr.write("Thank you! Couldn't reach the survey server just now; "
                         "circtools will retry silently on a later run.\n\n")
    _save_state(state)

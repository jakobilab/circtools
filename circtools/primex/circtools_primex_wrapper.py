#!/usr/bin/env python3


import sys
import primer3




def _na_row(seq_id: str, index: int) -> str:
    """Return a single NA row matching the original R output format."""
    label = f"{seq_id}_{index}"
    return "\t".join([label] + ["NA"] * 9)


def design_primers(
    seq_id: str,
    exon1_seq: str,
    exon2_seq: str,
    product_min: int,
    product_max: int,
    junction_mode: str,
    num_pairs: int,
) -> list[dict]:


    template = exon2_seq + exon1_seq
    bsj_pos  = len(exon2_seq)          # 0-based index of the BSJ


    if junction_mode == "f":
        # left_start, left_len, right_start, right_len
        # right_start / right_len = -1 means unconstrained in primer3-python
        ok_regions = [[bsj_pos - 15, 30, -1, -1]]
    elif junction_mode == "r":
        ok_regions = [[-1, -1, bsj_pos - 15, 30]]
    else:
        # left primer in exon2 (positions 0 .. bsj_pos-1)
        # right primer in exon1 (positions bsj_pos .. end)
        ok_regions = [[0, bsj_pos, bsj_pos, len(exon1_seq) - 1]]

    seq_args = {
        "SEQUENCE_ID":                         seq_id,
        "SEQUENCE_TEMPLATE":                   template,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": ok_regions,
    }

    global_args = {
        "PRIMER_TASK":                "generic",
        "PRIMER_PICK_LEFT_PRIMER":    1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_PICK_RIGHT_PRIMER":   1,
        "PRIMER_NUM_RETURN":          int(num_pairs),
        "PRIMER_PRODUCT_SIZE_RANGE":  [[product_min, product_max]],
    }

    try:
        result = primer3.design_primers(seq_args, global_args)
    except Exception as exc:
        sys.stderr.write(f"primer3 error for {seq_id}: {exc}\n")
        return []

    num_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if num_returned == 0:
        return []

    pairs = []
    for i in range(num_returned):
        left_pos  = result.get(f"PRIMER_LEFT_{i}",  [None, None])
        right_pos = result.get(f"PRIMER_RIGHT_{i}", [None, None])

        pairs.append({
            "PRIMER_LEFT_SEQUENCE":    result.get(f"PRIMER_LEFT_{i}_SEQUENCE",    "NA"),
            "PRIMER_RIGHT_SEQUENCE":   result.get(f"PRIMER_RIGHT_{i}_SEQUENCE",   "NA"),
            "PRIMER_LEFT":             f"{left_pos[0]},{left_pos[1]}"  if left_pos[0]  is not None else "NA",
            "PRIMER_RIGHT":            f"{right_pos[0]},{right_pos[1]}" if right_pos[0] is not None else "NA",
            "PRIMER_LEFT_TM":          result.get(f"PRIMER_LEFT_{i}_TM",          "NA"),
            "PRIMER_RIGHT_TM":         result.get(f"PRIMER_RIGHT_{i}_TM",         "NA"),
            "PRIMER_LEFT_GC_PERCENT":  result.get(f"PRIMER_LEFT_{i}_GC_PERCENT",  "NA"),
            "PRIMER_RIGHT_GC_PERCENT": result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", "NA"),
            "PRIMER_PAIR_PRODUCT_SIZE": result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", "NA"),
        })

    return pairs


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 5:
        sys.stderr.write(
            "Usage: circtools_primex_wrapper.py <data_file> "
            "<min_size,max_size> <junction_mode> <num_pairs>\n"
        )
        sys.exit(1)

    data_file_name = sys.argv[1]
    product_size   = [int(x) for x in sys.argv[2].split(",")]
    junction_mode  = sys.argv[3]
    num_pairs      = int(sys.argv[4])

    product_min, product_max = product_size[0], product_size[1]

    output_rows: list[str] = []

    with open(data_file_name, "r") as fh:
        for raw_line in fh:
            raw_line = raw_line.rstrip("\n")
            if not raw_line:
                continue

            columns = raw_line.split("\t")
            seq_id   = columns[0]
            exon2    = columns[1]   
            exon1    = columns[2] if len(columns) > 2 else ""  
            if not exon1:
                mid   = len(exon2) // 2
                exon1 = exon2[mid:]  
                exon2 = exon2[:mid]  

            pairs = design_primers(
                seq_id, exon1, exon2,
                product_min, product_max,
                junction_mode, num_pairs,
            )

            if pairs:
                for idx, pair in enumerate(pairs, start=1):
                    label = f"{seq_id}_{idx}"
                    row = "\t".join([
                        label,
                        str(pair["PRIMER_LEFT_SEQUENCE"]),
                        str(pair["PRIMER_RIGHT_SEQUENCE"]),
                        str(pair["PRIMER_LEFT"]),
                        str(pair["PRIMER_RIGHT"]),
                        str(pair["PRIMER_LEFT_TM"]),
                        str(pair["PRIMER_RIGHT_TM"]),
                        str(pair["PRIMER_LEFT_GC_PERCENT"]),
                        str(pair["PRIMER_RIGHT_GC_PERCENT"]),
                        str(pair["PRIMER_PAIR_PRODUCT_SIZE"]),
                    ])
                    output_rows.append(row)
            else:
                output_rows.append(_na_row(seq_id, 1))

    print("\n".join(output_rows))


if __name__ == "__main__":
    main()
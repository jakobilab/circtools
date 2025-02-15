.. circtools documentation master file, created by
   sphinx-quickstart on Thu Feb 22 15:40:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================================================================
circtools: a one-stop software solution for circular RNA research
==================================================================

Introduction
============

Circular RNAs (circRNAs) originate through back-splicing events from linear primary transcripts, are resistant to exonucleases, typically not polyadenylated, and have been shown to be highly specific for cell type and developmental stage. Although few circular RNA molecules have been shown to exhibit miRNA sponge function, for the vast majority of circRNAs however, their function is yet to be determined.

The prediction of circular RNAs is a multi-stage bioinformatics process starting with raw sequencing data and usually ending with a list of potential circRNA candidates which, depending on tissue and condition may contain hundreds to thousands of potential circRNAs. While there already exist a number of tools for the prediction process (e.g. DCC and circTest), publicly available downstream analysis tools are rare.

We developed **circtools**, a modular, Python3-based framework for circRNA-related tools that unifies several functionalities in single command line driven software. The command line follows the `circtools subcommand` standard that is employed in samtools or bedtools. Circtools includes modules for RBP enrichment screenings, circRNA primer design, as well as interfaces to the processing tools FUCHS and DCC; miRNA seed analysis and differential exon usage will we available in the upcoming release.

We intend to add more and more modules in the future in order to provide a comprehensive bioinformatics toolbox and also encourage users to contribute modules to circtools.

**circtools** is developed at the `Jakobi Lab <https://jakobilab.org>`_, part of the `Translational Cardiovascular Research Center (TCRC) <https://https://phoenixmed.arizona.edu/tcrc/>`_ in the Department of Internal Medicine at `The University of Arizona College of Medicine – Phoenix <https://https://phoenixmed.arizona.edu/>`_ with contributions from the `Dieterich Lab <https://dieterichlab.org>`_ at the `University Hospital Heidelberg <https://www.heidelberg-university-hospital.com/>`_, where `Tobias Jakobi <https://tobias.jako.bi>`_ initially developed the software.

Citation
==========

If you are using **circtools** in your work, **please cite the following publication** to help us provide support and continued development for the software:

`Jakobi, T., Uvarovskii, A., & Dieterich, C. (2019). circtools — a one-stop software solution for circular RNA research. Bioinformatics, 35(13), 2326-2328. <https://academic.oup.com/bioinformatics/article/35/13/2326/5194340>`_

License
========
**circtools** is freely available under the GNU General Public License v3.0.

Issues
===============

Problems or issues should be reported `via the GitHub issue system <https://github.com/jakobilab/circtools/issues>`_.

Table of contents
==================

.. toctree::
   :maxdepth: 3

   Installation.rst
   Detect.rst
   Nanopore.rst
   Reconstruct.rst
   Conservation.rst
   Padlock.rst
   Quickcheck.rst
   Circtest.rst
   Enrichment.rst
   Exon.rst
   Primer.rst
   Support.rst
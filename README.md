# Proteomics Research Project

## Motivation:
Traditionally, protein abundance levels are used to detect intermolecular changes and interactions, but this method is insufficient for dealing with chemical modifications that do not affect protein levels. This Proteomics Data Science (protds) project tests an alternative approach of using protein structural states to predict modifications.

## About:
This repository contains a collection of code used to analyze protein data obtained from mass spectrometry. The main component of the project is a program named protds that collects and organizes protein structure information from multiple databases to be used for visualization and data analysis. These databases include:  
- Protein Data Bank (https://www.rcsb.org): source for PDB structure data
- MMTF (https://mmtf.rcsb.org): source for efficient protein structures in binary format
- UniProtKB (https://www.uniprot.org): source for protein binding sites and feature data
- AlphaFoldDB (https://alphafold.ebi.ac.uk): source for predicted protein structures

From an input dataset of proteins, the program pulls information for each entry and saves it as a dictionary that can be cached as a pickle file for future access. Examples of what can be done with this information is demonstrated in the Jupyter Notebook files provided:
- protds_v3.ipynb compares and highlights modified peptide sequences and known binding sites
- sequences.ipynb applies regression analysis to generate planes through a set of peptides 

Python scripts for other proteomics-related projects are stored in this repository as well. These include code for data preprocessing, GRAVY index calculations, and PDB result counting.
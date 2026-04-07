#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Erin Nicole Decocker
# edecocke@charlotte.edu
# ID: 801442694

"""
merge_practice.py

Merges gene expression and mutation datasets, computes mean expression,
and performs correlation analysis to test whether highly expressed genes
tend to accumulate more mutations.
"""

import pandas as pd


# ---------------------------------------------------------------------------
# Classes
# ---------------------------------------------------------------------------

class TheMutantMashup:
    """
    Loads, merges, and analyzes gene expression and mutation datasets.

    Merging is performed on the shared 'Gene' column using pd.merge(),
    which aligns rows by a specified key column rather than by index.
    This is the appropriate method when the shared key is a regular
    column rather than the DataFrame index.
    """

    def __init__(self, expression_data: dict, mutation_data: dict):
        """
        Initialises both input DataFrames from raw dictionaries.

        Parameters
        ----------
        expression_data : dict
            Dictionary containing gene expression values across tissues.
            Must include a 'Gene' key for merging.
        mutation_data : dict
            Dictionary containing mutation counts per gene.
            Must include a 'Gene' key for merging.

        Ensures
        -------
        self.df_expression and self.df_mutations are valid pandas
        DataFrames ready for merging. self.df_merged is initialised
        as None until merge_datasets() is called.

        Returns
        -------
        None
        """
        self.df_expression = pd.DataFrame(expression_data)
        self.df_mutations = pd.DataFrame(mutation_data)
        self.df_merged = None

    def merge_datasets(self) -> pd.DataFrame:
        """
        Merges the expression and mutation DataFrames on the 'Gene' column
        using an inner join and prints the result.

        pd.merge() is used here rather than .join() because the shared
        key ('Gene') is a regular column, not the DataFrame index.
        An inner join keeps only genes present in both tables.

        Parameters
        ----------
        None

        Ensures
        -------
        self.df_merged is assigned the merged DataFrame.
        Only genes present in both input DataFrames are retained.

        Returns
        -------
        pd.DataFrame
            The merged DataFrame containing expression and mutation columns.
        """
        self.df_merged = pd.merge(self.df_expression, self.df_mutations, on='Gene')
        print(f"\n{'-'*60}")
        print("Merged DataFrame")
        print(f"{'-'*60}")
        print(self.df_merged)
        return self.df_merged

    def compute_mean_expression(self) -> None:
        """
        Computes mean expression across all three tissue columns for each
        gene and appends the result as a new 'Mean_Expression' column.

        Averaging is performed across columns (axis=1) so that each gene
        receives a single mean value representing its overall expression level.

        Parameters
        ----------
        None

        Ensures
        -------
        self.df_merged is modified in place to include a 'Mean_Expression'
        column. merge_datasets() must be called before this method.

        Returns
        -------
        None
        """
        tissue_cols = ['Expression_Heart', 'Expression_Liver', 'Expression_Brain']
        self.df_merged['Mean_Expression'] = self.df_merged[tissue_cols].mean(axis=1).round(3)
        print(f"\n{'-'*60}")
        print("DataFrame with Mean Expression")
        print(f"{'-'*60}")
        print(self.df_merged)

    def unleash_the_correlation_matrix(self) -> pd.DataFrame:
        """
        Computes the Pearson correlation matrix for all numeric columns
        including tissue expression values, mean expression, and mutation counts.

        DataFrame.corr() computes pairwise Pearson correlation coefficients.
        Values range from -1 (perfect negative correlation) to +1 (perfect
        positive correlation). A value close to +1 between Mean_Expression
        and Mutations would suggest that highly expressed genes tend to
        accumulate more mutations.

        Parameters
        ----------
        None

        Ensures
        -------
        The correlation matrix includes all numeric columns in self.df_merged.
        compute_mean_expression() must be called before this method.

        Returns
        -------
        pd.DataFrame
            A square correlation matrix with variables as both row and column labels.
        """
        cols = [
            'Expression_Heart',
            'Expression_Liver',
            'Expression_Brain',
            'Mean_Expression',
            'Mutations'
        ]
        correlation_matrix = self.df_merged[cols].corr()
        print(f"\n{'-'*60}")
        print("Correlation Matrix")
        print(f"{'-'*60}")
        print(correlation_matrix.round(3))
        return correlation_matrix

    def export_to_csv(self, filename: str = "gene_expression_mutations.csv") -> None:
        """
        Exports the merged and annotated DataFrame to a CSV file.

        Parameters
        ----------
        filename : str
            Output file path. Defaults to 'gene_expression_mutations.csv'
            in the current working directory.

        Ensures
        -------
        The CSV file is written without the integer row index, as the
        'Gene' column already serves as a meaningful identifier.
        merge_datasets() and compute_mean_expression() must be called
        before this method.

        Returns
        -------
        None
        """
        # index=False omits the default integer row index from the output file
        self.df_merged.to_csv(filename, index=False)
        print(f"\n Results saved to '{filename}'")


# ---------------------------------------------------------------------------
# Reflection
# ---------------------------------------------------------------------------

def interpret_the_chaos() -> None:
    """
    Prints an interpretation of the correlation results
    """
    interpretation = (
        "A high positive correlation (close to +1) between Mean_Expression\n"
        "   and Mutations suggests that genes with higher expression levels\n"
        "   also tend to accumulate more mutations. This pattern can occur\n"
        "   in actively transcribed genomic regions, where increased\n"
        "   replication and transcription activity raises the likelihood\n"
        "   of replication errors or DNA damage.\n\n"
        "   Note: correlation does not imply causation. Additional\n"
        "   statistical testing would be required before drawing\n"
        "   biological conclusions from this data."
    )

    print(f"\n{'-'*60}")
    print("Interpretation of Results")
    print(f"{'-'*60}")
    print(f"\n   {interpretation}")

if __name__ == "__main__":

    expression_data = {
        'Gene': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
        'Expression_Heart': [5.2, 3.3, 7.1, 4.5, 2.8],
        'Expression_Liver': [4.8, 2.9, 6.5, 4.1, 2.5],
        'Expression_Brain': [6.1, 3.7, 8.2, 5.3, 3.0]
    }

    mutation_data = {
        'Gene': ['BRCA1', 'TP53', 'EGFR', 'MYC', 'KRAS'],
        'Mutations': [12, 8, 15, 10, 6]
    }
    mashup = TheMutantMashup(
        expression_data=expression_data,
        mutation_data=mutation_data
    )

    mashup.merge_datasets()
    mashup.compute_mean_expression()
    mashup.unleash_the_correlation_matrix()
    interpret_the_chaos()
    mashup.export_to_csv("gene_expression_mutations.csv")

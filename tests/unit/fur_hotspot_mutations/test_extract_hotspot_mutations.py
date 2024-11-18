import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

from fur_hotspot_mutations.extract_hotspot_mutations import extract_hotspot_mutations

# Sample mock data for testing
mock_maf_data = pd.DataFrame(
    {
        "Hugo_Symbol": ["TP53", "TP53", "TP53", "BRCA1", "BRCA1", "BRCA1", "EGFR"],
        "Chromosome": ["17", "17", "17", "17", "17", "17", "7"],
        "Start_Position": [
            7574000,
            7574000,
            7574000,
            43045629,
            43045629,
            43045629,
            55249071,
        ],
        "End_Position": [
            7574001,
            7574001,
            7574001,
            43045630,
            43045630,
            43045630,
            55249072,
        ],
        "Reference_Allele": ["C", "C", "C", "A", "A", "A", "G"],
        "Tumor_Seq_Allele2": ["T", "T", "T", "T", "T", "T", "A"],
        "Tumor_Sample_Barcode": [
            "Sample1",
            "Sample2",
            "Sample3",
            "Sample1",
            "Sample2",
            "Sample3",
            "Sample1",
        ],
    }
)

# Expected output (EGFR mutation is filtered out because it occurs in only one sample)
expected_output_data = pd.DataFrame(
    {
        "Hugo_Symbol": ["TP53", "TP53", "TP53", "BRCA1", "BRCA1", "BRCA1"],
        "Chromosome": ["17", "17", "17", "17", "17", "17"],
        "Start_Position": [7574000, 7574000, 7574000, 43045629, 43045629, 43045629],
        "End_Position": [7574001, 7574001, 7574001, 43045630, 43045630, 43045630],
        "Reference_Allele": ["C", "C", "C", "A", "A", "A"],
        "Tumor_Seq_Allele2": ["T", "T", "T", "T", "T", "T"],
        "Tumor_Sample_Barcode": [
            "Sample1",
            "Sample2",
            "Sample3",
            "Sample1",
            "Sample2",
            "Sample3",
        ],
    }
)


# Mock pd.read_csv and to_csv
@patch("pandas.read_csv")
@patch("pandas.DataFrame.to_csv")
def test_extract_hotspot_mutations(mock_to_csv, mock_read_csv):
    # Setup the mock for pd.read_csv
    mock_read_csv.return_value = mock_maf_data

    # Call the function
    output = extract_hotspot_mutations("input_maf_file.maf", "output_maf_file.maf")

    # Assert pd.read_csv was called correctly
    mock_read_csv.assert_called_once_with(
        "input_maf_file.maf", sep="\t", comment="#", low_memory=False
    )

    # Assert to_csv was called with correct arguments
    mock_to_csv.assert_called_once_with("output_maf_file.maf", sep="\t", index=False)

    # Verify the output DataFrame matches the expected result
    pd.testing.assert_frame_equal(output, expected_output_data)


def test_missing_required_columns():
    # Create mock data with missing columns
    mock_data = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53"],  # Missing other required columns
            "Chromosome": ["17"],
        }
    )

    with patch("pandas.read_csv", return_value=mock_data):
        with pytest.raises(
            ValueError,
            match="Required column 'Start_Position' not found in the MAF file.",
        ):
            extract_hotspot_mutations("input_maf_file.maf", "output_maf_file.maf")


def test_no_hotspot_mutations():
    # Create mock data where no mutation occurs in >2 samples
    mock_data = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "BRCA1"],
            "Chromosome": ["17", "17"],
            "Start_Position": [7574000, 43045629],
            "End_Position": [7574001, 43045630],
            "Reference_Allele": ["C", "A"],
            "Tumor_Seq_Allele2": ["T", "T"],
            "Tumor_Sample_Barcode": [
                "Sample1",
                "Sample1",
            ],  # Each mutation only in 1 sample
        }
    )

    with patch("pandas.read_csv", return_value=mock_data):
        output = extract_hotspot_mutations(
            "input_maf_file.maf", "output_maf_file.maf", save_to_file=False
        )

        # Output should be an empty DataFrame
        assert output.empty


def test_duplicate_sample_barcodes():
    # Create mock data with additional mutations
    mock_data = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "TP53", "TP53", "EGFR", "EGFR", "BRCA1", "BRCA1"],
            "Chromosome": ["17", "17", "17", "7", "7", "17", "17"],
            "Start_Position": [
                7574000,
                7574000,
                7574000,
                550000,
                550000,
                430000,
                430000,
            ],
            "End_Position": [7574001, 7574001, 7574001, 550001, 550001, 430001, 430001],
            "Reference_Allele": ["C", "C", "C", "A", "A", "G", "G"],
            "Tumor_Seq_Allele2": ["T", "T", "T", "G", "G", "A", "A"],
            "Tumor_Sample_Barcode": [
                "Sample1",
                "Sample1",
                "Sample2",  # TP53: meets criteria
                "Sample3",
                "Sample4",  # EGFR: meets criteria
                "Sample5",
                "Sample5",  # BRCA1: does NOT meet criteria
            ],
        }
    )

    with patch("pandas.read_csv", return_value=mock_data):
        output = extract_hotspot_mutations(
            "input_maf_file.maf", "output_maf_file.maf", save_to_file=False
        )

        # Check that the output is not empty
        assert not output.empty

        # Check the number of unique mutations that meet the criteria
        expected_hotspot_genes = ["TP53", "EGFR"]  # BRCA1 should not meet the criteria
        assert set(output["Hugo_Symbol"].unique()) == set(expected_hotspot_genes)

        # Ensure correct counts for individual genes
        tp53_samples = output[output["Hugo_Symbol"] == "TP53"]["Tumor_Sample_Barcode"]
        assert tp53_samples.nunique() == 2  # TP53 has 2 unique samples

        egfr_samples = output[output["Hugo_Symbol"] == "EGFR"]["Tumor_Sample_Barcode"]
        assert egfr_samples.nunique() == 2  # EGFR has 2 unique samples

        # Ensure BRCA1 does not appear in the output
        assert "BRCA1" not in output["Hugo_Symbol"].unique()

        # Ensure all rows in the output correspond to valid hotspots (deduplicated)
        assert len(output) == 4  # 2 for TP53 (deduplicated), 2 for EGFR

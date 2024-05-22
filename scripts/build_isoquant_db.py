import sys
import gffutils

def create_gff_database(gtf_path, db_path):
    """
    Creates a gffutils database from a GTF file.

    Args:
        gtf_path (str): Path to the input GTF file.
        db_path (str): Path to the output database file.
    """
    try:
        gffutils.create_db(
            gtf_path, 
            db_path, 
            force=True,        # Overwrite existing database if it exists
            keep_order=True,   # Maintain the order of features from the GTF
            merge_strategy='error', # Raise an error if features overlap
            sort_attribute_values=True, # Sort attribute values for consistency
            disable_infer_transcripts=True,  
            disable_infer_genes=True,       
        )
        print(f"Database created successfully at: {db_path}")
    except Exception as e:
        print(f"Error creating database: {e}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <gtf_file> <output_db_file>")
        sys.exit(1)  # Exit with an error code

    gtf_path = sys.argv[1]
    db_path = sys.argv[2]

    create_gff_database(gtf_path, db_path)
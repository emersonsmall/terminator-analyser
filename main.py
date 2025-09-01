import gffutils

def main():
    gffutils.create_db("../genomes/thaliana.gff", dbfn="thaliana_gff.db", keep_order=True, merge_strategy="create_unique", sort_attribute_values=True)
    print("Database created.")
    return

if __name__ == "__main__":
    main()

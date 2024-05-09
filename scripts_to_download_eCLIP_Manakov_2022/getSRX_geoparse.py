# Importing libraries
import sys
import GEOparse

# Retrieving relevant SRXs from GSE in GEO DB
def generate_srx_list(geo_series):
    gse = GEOparse.get_GEO(geo=geo_series, destdir="./")
    srx_list = []
    substrings_to_check = ["MusLiver", "Enriched", "oe"]
    for gsm_name, gsm in gse.gsms.items():
        title = gsm.metadata['title'][0]
        if "ChimeCLIP" in title:
            if not any(substring in title for substring in substrings_to_check):
                sra_link = gsm.metadata['relation'][1]
                srx_list.append(
                    sra_link[sra_link.index('=') + 1:]
                )                                  
    space_sep_srx_list = " ".join(srx_list)
    return space_sep_srx_list

if __name__ == "__main__":
    geo_id = sys.argv[1]
    print(generate_srx_list(geo_id)) 

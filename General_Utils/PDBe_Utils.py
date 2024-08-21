from tempfile import gettempdir
from typing import List, Dict, Tuple

import biotite.database.rcsb as rcsb
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import numpy as np
import pandas as pd


class PDBe_Search:
    BASE_URL = "https://www.ebi.ac.uk/pdbe/"
    SEARCH_URL = BASE_URL + 'search/pdb/select?'

    @classmethod
    def make_request_post(cls, search_dict, number_of_rows):
        """
        Makes a POST request to the PDBe API
        """
        import requests

        if 'rows' not in search_dict:
            search_dict['rows'] = number_of_rows

        search_dict['wt'] = 'json'
        response = requests.post(cls.SEARCH_URL, data=search_dict)

        if response.status_code == 200:
            return response.json()

        else:
            print(f"[No data retrieved - {response.status_code}] {response.text}")

        return {}

    @classmethod
    def format_search_terms_post(cls, search_terms, filter_terms=None, **kwargs):
        """
        Formats the search terms for the PDBe API
        """
        # Variable to return
        return_variables = {'q': str(search_terms)}

        if filter_terms:
            fl = ','.join(filter_terms)
            return_variables['fl'] = fl

        for arg in kwargs:
            return_variables[arg] = kwargs[arg]

        return return_variables

    @classmethod
    def run_search(cls, search_terms, filter_terms=None, number_of_rows=500, **kwargs):
        """
        Run the search with a set of search terms
        """
        search_params = cls.format_search_terms_post(search_terms=search_terms, filter_terms=filter_terms, **kwargs)

        if search_params:
            response = cls.make_request_post(search_dict=search_params, number_of_rows=number_of_rows)

            if response:
                results = response.get('response', {}).get('docs', [])
                print(f'Number of results for {search_terms}: {len(results)}')
                return results

        print('No results')
        return []

    @classmethod
    def pandas_dataset(cls, list_of_results):
        """
        Updates lists to strings for loading into Pandas
        """
        for row in list_of_results:
            for data in row:
                if type(row[data]) == list:
                    # If there are any numbers in the list change them into strings
                    row[data] = [str(a) for a in row[data]]

                    # Unique and sort the list and then change the list into a string
                    row[data] = ','.join(sorted(list(set(row[data]))))

        df = pd.DataFrame(list_of_results)
        return df


# # Example Usage
# search_terms = Q(uniprot_accession='P01116', mutation="y")
# filter_terms = ['pdb_id', "title", "mutation"]
#
# # Run the search
# pdbe_result = PDBe_Search.run_search(search_terms, filter_terms)
# pdbe_df = PDBe_Search.pandas_dataset(pdbe_result)

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

def pdbe_infos(pdbid, case_search_url: str) -> Dict:
    import requests
    BASE_URL = "https://www.ebi.ac.uk/pdbe/"
    SEARCH_URL = BASE_URL + case_search_url

    # Ask to the PDBe
    response = requests.get(SEARCH_URL + pdbid).json()
    infos = response[pdbid]
    return infos

def pdbe_spot_title(pdbid) -> Dict:
    import requests
    BASE_URL = "https://www.ebi.ac.uk/pdbe/"
    SEARCH_URL = BASE_URL + 'api/pdb/entry/summary/'

    # Ask to the PDBe
    response = requests.get(SEARCH_URL + pdbid).json()
    infos = response[pdbid]
    return infos[0]


def pdbe_spot_resids(pdbid) -> List:
    import requests
    BASE_URL = "https://www.ebi.ac.uk/pdbe/"
    SEARCH_URL = BASE_URL + 'api/pdb/entry/ligand_monomers/'

    # Ask to the PDBe
    response = requests.get(SEARCH_URL + pdbid).json()
    infos = response[pdbid]

    resids: List = [item["chem_comp_id"] for item in infos]
    excluded_res_list = ["GTP", "ARG", "SPM", "NCO", "SO4", "PO4", "CL",
                         "3AD", "IPA", "1PE", "IRI", "ACT", "SPD", "GOL",
                         "NAD", "PG4", "SPM", "DMS", "SPK", "FRU", "MG",
                         "GLC", "SIN", "CAC", "NH4", "G3A", "EDO", "GDP", "CA"]

    resids: List = list(np.unique([res for res in resids if res not in excluded_res_list]))

    return resids


def biotite_spot_resids(pdbid) -> List:
    protein_path = rcsb.fetch(pdbid, "pdb", gettempdir(), overwrite=True)
    reader = pdb.PDBFile.read(protein_path)
    struct_array = reader.get_structure(model=1)

    crystal = struct_array[
        ~(struc.filter_amino_acids(struct_array) |
          struc.filter_nucleotides(struct_array) |
          struc.filter_monoatomic_ions(struct_array) |
          struc.filter_solvent(struct_array))
    ]
    excluded_res_list = ["GTP", "ARG", "SPM", "NCO", "SO4", "PO4", "CL",
                         "3AD", "IPA", "1PE", "IRI", "ACT", "SPD", "GOL",
                         "NAD", "PG4", "SPM", "DMS", "SPK", "FRU", "MG",
                         "GLC", "SIN", "CAC", "NH4", "G3A", "EDO", "GDP", "CA"]
    crystal: List = list(np.unique(crystal[~np.isin(crystal.res_name, excluded_res_list)].res_name))

    return crystal


def pdbe_spot_smile(resid) -> str:
    import requests
    BASE_URL = "https://www.ebi.ac.uk/pdbe/"
    SEARCH_URL = BASE_URL + 'api/pdb/compound/summary/'

    # Ask to the PDBe
    response = requests.get(SEARCH_URL + resid).json()
    infos = response[resid]

    smiles = infos[0]["smiles"][-1]["name"]
    return smiles

import sys
import os
import pickle
from deepARG_psilibs.features_vectorizer import Features_extractor
from deepARG_psilibs.create_directories import create_dir
from deepARG_psilibs.network import train_network
from deepARG_psilibs.prokka import header_reducer, prokka_translator

class Controller(object):

    """Manage the actions of the subcommands.
    The Controller takes care of providing the arguments like path
    names and the parallel processing of tasks.
    """

    def __init__(self, args):
        """Create an instance."""
        self._args = args
    
    def create_project_folder(self):
        create_dir(self._args.root_path)

    def vectorize_features(self):
        ft = Features_extractor(self._args.fasta_path,
                                self._args.root_path)
        data = ft.PSSM_comp_Extractor()
        pickle.dump(data,
                    open("{}/input/{}.pkl".format(self._args.root_path,
                                                  self._args.dataset_name),
                         "wb")
                    )

    def create_model(self):
        dataset = pickle.load(open(self._args.dataset_path,
                              "rb"))
        netx = train_network(dataset[0])
        netx.fit(dataset[0], dataset[1])
        filename = "{}/model/{}.pkl".format(self._args.root_path,
                                            self._args.model_name)
        pickle.dump(netx, open(filename, 'wb'))

    def predict_ABG(self):
        if self._args.prot_type is True:
            protein_seq = self._args.fasta_sequence
        else:
            shrt_hdr = header_reducer(self._args.root_path,
                                  self._args.fasta_sequence,
                                  self._args.outfile_prefix)
            protein_seq = prokka_translator(shrt_hdr,
                                        self._args.root_path,
                                        self._args.outfile_prefix)
            if not os.path.exists(protein_seq):
                sys.stdout.write("Prokka did not create protein sequence \nEXITING!!! \n")
                sys.exit(2)

        ft = Features_extractor(protein_seq,
                               self._args.root_path)
        ft.PSSM_comp_Extractor_predict(self._args.model_path,
                                       "{}/predictions/{}.csv".format(
                                                        self._args.root_path,
                                                        self._args.outfile_prefix))





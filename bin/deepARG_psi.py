#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
import argparse
from deepARG_psilibs.controller import Controller


def main():
    parser = create_parser()
    args = parser.parse_args()
    controller = Controller(args)
    args.func(controller)

def create_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="commands")
    
    # create directories
    create_dir_parser = subparsers.add_parser(
        "create", help="Create root directory")
    create_dir_parser.add_argument(
        "--root_path", "-rp", help="path of root folder")
    create_dir_parser.set_defaults(func=create_dir)

    # feature vectorizer
    feature_vectorizer_parser = subparsers.add_parser(
        "vectorize", help="Vectorize a set of sequences")
    feature_vectorizer_parser.add_argument(
        "--fasta_path", "-fp", help="Path to fasta file")
    feature_vectorizer_parser.add_argument(
        "--root_path", "-rp", help="Path to root directory")
    feature_vectorizer_parser.add_argument(
        "--dataset_name", "-dn", help="prexix to name dataset")
    feature_vectorizer_parser.set_defaults(func=vectorize_features)
    
    # Train model
    train_model_parser = subparsers.add_parser(
        "train", help="train model on dataset")
    train_model_parser.add_argument(
        "--dataset_path", "-dp", help="Path to fasta file")
    train_model_parser.add_argument(
        "--model_name", "-mn", help="Name to call model")
    train_model_parser.add_argument(
        "--root_path", "-rp", help="Path to root directory")
    train_model_parser.set_defaults(func=create_model)

    # predictor caller
    prediction_parser = subparsers.add_parser(
        "predict", help="predict antibiotic category")
    prediction_parser.add_argument(
        "--prediction_type", "-pet", help="Prediction type CTD or PSSM")
    prediction_parser.add_argument(
        "--fasta_sequence", "-fs", help="Path to fasta file")
    prediction_parser.add_argument(
        "--prot_type", "-pt", default=False,
        action="store_true", help="Specify protein sequence type")
    prediction_parser.add_argument(
        "--model_path", "-mp", help="Path to neural net model")
    prediction_parser.add_argument(
        "--outfile_prefix", "-op", help="Prefix for out file")
    prediction_parser.add_argument(
        "--root_path", "-rp", help="path of root folder")
    prediction_parser.set_defaults(func=predict_ABG)
    return parser


def create_dir(controller):
    controller.create_project_folder()


def vectorize_features(controller):
    controller.vectorize_features()


def create_model(controller):
    controller.create_model()


def predict_ABG(controller):
    controller.predict_ABG()

if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
from lasagne import layers
from lasagne.updates import sgd, nesterov_momentum
from nolearn.lasagne import NeuralNet
import numpy as np
import theano
from deepARG_psilibs.model import flex_model


def train_network(xtrain):
    network = NeuralNet(
        layers=flex_model(xtrain.shape[1], 400, 45),
        update=nesterov_momentum,
        update_learning_rate=0.01,
        update_momentum=0.9,
        regression=False,
        max_epochs=40,
        eval_size=0.2,
        verbose=1,
    )

    return network

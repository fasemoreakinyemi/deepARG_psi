#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
from lasagne import layers
import lasagne

def model(input_size, output_size):
    # input size: x_train.shape[1]
    # output size: len(set(Y)
    return [
            (layers.InputLayer,  {'shape':(None, input_size)}),
            
            (layers.DenseLayer,  {'num_units':148}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':74}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':37}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':20}),
            
            (layers.DenseLayer, {'num_units':output_size,
                                    'nonlinearity':lasagne.nonlinearities.softmax
                                }),
        ]

def flex_model(input_size, unit_size, output_size):
    # input size: x_train.shape[1]
    # output size: len(set(Y)
    return [
            (layers.InputLayer,  {'shape':(None, input_size)}),
            
            (layers.DenseLayer,  {'num_units':unit_size}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':int(unit_size*0.5)}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':int(unit_size*0.25)}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':int(unit_size*0.05)}),
            
            (layers.DenseLayer, {'num_units':output_size,
                                    'nonlinearity':lasagne.nonlinearities.softmax
                                }),
        ]

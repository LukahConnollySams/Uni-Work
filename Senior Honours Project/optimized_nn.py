import tensorflow as tf
import numpy as np

def optimal_layers(max_layers, x_train, y_train, x_valid, y_valid, node_size=100, activation='relu', batch_size=16, epochs=2000, patience=15):
    """
    Creates multiple neural networks with different numbers of hidden layers. Then trains and evaluates them, returning the loss value for each hidden layer number.

    Parameters
    ----------
    max_layers : int
        Max number of hidden layers to train
    x_train : array-like
        Input training data
    y_train : array-like
        Output training data
    x_valid : array-like
        Input validation data
    y_valid : array-like
        Output validation data
    node_size : int, optional
        Number of nodes per hidden layer, by default 100
    activation : str, optional
        Activation function for hidden layers, by default 'relu'
    batch_size : int, optional
        Batch size for evaluation, by default 16
    epochs : int, optional
        Number of epochs for training, by default 2000
    patience : int, optional
        Patience for EarlyStopping callback function, by default 15

    Returns
    -------
    array-like/list
        list of layers values and corresponding loss values.
    """

    #arrays/list for final values
    layer_number = np.arange(1, max_layers + 1)
    layers_losses = []
    
    #loop for the number of layers
    for i in range(len(layer_number)):
        
        #initialise model
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Input(shape=len(x_train[0])))

        #loop to add set amount of layers
        x = 0
        while x < layer_number[i]:

            model.add(tf.keras.layers.Dense(node_size, activation=activation))
            x += 1

        #add output layer
        model.add(tf.keras.layers.Dense(len(y_train[0]), activation='linear'))
        
        #compile model with adam optimizer and MeanAbsoluteError loss function
        model.compile(optimizer=tf.keras.optimizers.Adam(), loss=tf.keras.losses.MeanAbsoluteError(), metrics=['accuracy'])
        
        #train model with Early Stopping callback
        model.fit(x_train, y_train, epochs=epochs, validation_data=(x_valid, y_valid), callbacks=[tf.keras.callbacks.EarlyStopping(monitor='loss', patience=patience)])

        #evaluate model and return history
        history = model.evaluate(x_valid, y_valid, batch_size=batch_size, return_dict=True)

        #record final loss value for model
        layers_losses.append(history['loss'])

        #clear for new model
        tf.keras.backend.clear_session()


    return layer_number, layers_losses

def optimal_nodes(max_nodes, x_train, y_train, x_valid, y_valid, start_nodes=100, node_interval=100, layer_number=2, activation='relu', batch_size=16, epochs=2000, patience=15):
    """
    Creates multiple neural networks with different numbers of nodes per hidden layer. Then trains and evaluates them, returning the loss value for each node per hidden layer value.

    Parameters
    ----------
    max_nodes : int
        Maximum number of nodes per hidden layer
    x_train : array-like
        Input training data
    y_train : array-like
        Output training data
    x_valid : array-like
        Input validation data
    y_valid : array-like
        Output validation data
    start_nodes : int, optional
        Starting point for different model's nodes per hidden layer, by default 100
    node_interval : int, optional
        Nodes per hidden layer interval for different models, by default 100
    layer_number : int, optional
        Number of hidden layers, by default 2
    activation : str, optional
        Activation function for hidden layers, by default 'relu'
    batch_size : int, optional
        Batch size for evaluation, by default 16
    epochs : int, optional
        Number of epochs for training, by default 2000
    patience : int, optional
        Patience for EarlyStopping callback function, by default 15

    Returns
    -------
    array-like/list
        list of nodes per hidden layer values and corresponding loss values.
    """

    #arrays/list for final values
    node_number = np.arange(start_nodes, max_nodes + node_interval, node_interval)
    nodes_losses = []
    
    #loop for the number of nodes per hidden layer
    for i in range(len(node_number)):
        
        #initialise model
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Input(shape=(len(x_train[0]),)))
        
        #loop to add set amount of layers
        x = 0
        while x < layer_number:
            
            model.add(tf.keras.layers.Dense(node_number[i], activation=activation))
            x += 1

        #add output layer
        model.add(tf.keras.layers.Dense(len(y_train[0]), activation='linear'))

        #compile model with adam optimizer and MeanAbsoluteError loss function
        model.compile(optimizer=tf.keras.optimizers.Adam(), loss=tf.keras.losses.MeanAbsoluteError(), metrics=['accuracy'])
        
        #train model with Early Stopping callback
        model.fit(x_train, y_train, epochs=epochs, validation_data=(x_valid, y_valid), callbacks=[tf.keras.callbacks.EarlyStopping(monitor='loss', patience=patience)])

        #evaluate model and return history
        history = model.evaluate(x_valid, y_valid, batch_size=batch_size, return_dict=True)

        #record final loss value for model
        nodes_losses.append(history['loss'])

        #clear for new model
        tf.keras.backend.clear_session()
        
        
    return node_number, nodes_losses

def optimal_nn(nodes, layers, x_train, y_train, x_valid, y_valid, activation='relu', epochs=2000, patience=15):
    """
    Creates a neural network of specified size, then proceeding to train and test model from samples given, using Adam optimizer and MeanAbsolutError loss function.

    Parameters
    ----------
    nodes : int
        Number of nodes per hidden layer
    layers : int
        Number of hidden layers
    x_train : array-like
        Input training data
    y_train : array-like
        Output training data
    x_valid : array-like
        Input validation data
    y_valid : array-like
        Output validation data
    activation : str, optional
        Type of activation function used for hidden layers, by default 'relu'
    epochs : int, optional
        Number of epochs for training, by default 2000
    patience : int, optional
        Patience for EarlyStopping callback function, by default 15

    Returns
    -------
    tf.keras.Model
        Trained tensorflow model
    """

    #initialise model and input layer
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Input(shape=(len(x_train[0],))))
    
    #loop to add set amount of layers
    x = 0
    while x < layers:

        model.add(tf.keras.layers.Dense(nodes, activation=activation))
        x += 1
    
    #add output layer
    model.add(tf.keras.layers.Dense(len(y_train[0]), activation='linear'))
    
    #compile with Adam optimiser, Mean absolute error loss, and accuracy metrics
    model.compile(optimizer=tf.keras.optimizers.Adam(), loss=tf.keras.losses.MeanAbsoluteError(), metrics=['accuracy'])
    
    #train model
    model.fit(x_train, y_train, epochs=epochs, validation_data=(x_valid, y_valid), callbacks=[tf.keras.callbacks.EarlyStopping(monitor='loss', patience=patience)])

    #return model
    return model
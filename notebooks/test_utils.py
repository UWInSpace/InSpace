# This is a file of test functions for all of the functions in our class project.

#################################################################################
# Functions to run the UWInSpace tool for ChemE 545/546 in Winter 2021

def get_seq(email, prot_accession_num):
    from Bio import Entrez
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    
    Entrez.email = email

    gis = [prot_accession_num] 
    request = Entrez.epost("protein",id=",".join(map(str,gis)))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="protein",retmode="xml", webenv=webEnv, query_key=queryKey) 

    for r in Entrez.parse(handle):
        try:
            gi=int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
        except ValueError:
            gi=None
        #print(r['GBSeq_sequence'])
        seq = r['GBSeq_sequence']
    return seq



def usr_seq(user_input, email):
    """If user input is an integer (EID and not a sequence that would be a string), get the sequence
    """
    t = type(user_input)
    if t is int:
        user_sequence = get_seq(email, user_input)
    else:
        user_sequence = user_input
    
    return user_sequence



def get_prot_mass(AASeq):
    
    '''Tabulate approximate total mass of input amino acid string based on the sum of
    the average masses for every amino acid in the string.'''
    
    # Create dict with amino acid sybols as keys and associated average masses as values
    # Average masses taken from http://proteomicsresource.washington.edu/protocols06/masses.php
    molecular_mass_dict = {
    'A':71.0779, 
    'R':156.18568, 
    'G':57.05132, 
    'S':87.0773, 
    'P':97.11518, 
    'V':99.13106, 
    'T':101.10388, 
    'C':103.1429, 
    'L':113.15764, 
    'I':113.15764, 
    'N':114.10264, 
    'D':115.0874, 
    'Q':128.12922, 
    'K':128.17228, 
    'E':129.11398, 
    'M':131.19606, 
    'H':137.13928, 
    'F':147.17386, 
    'U':150.3079, 
    'Y':163.17326, 
    'W':186.2099, 
    'O':237.29816
    }
    
    # Calculate approx mass for amino acid sequence
    mass_list = []
    
    # Capitalize all letters in amino acid string prior to taking mass
    AASeq = AASeq.upper()
    
    # Raise error if amino acid sequence contains characters not in the dict
    acceptable_amino_acids_list = list(molecular_mass_dict.keys())

    for char in AASeq:
        if char not in acceptable_amino_acids_list:
            if char == ' ':
                continue
            else:
                raise ValueError('AASeq contains unacceptable characters')
    
    for amino_acid, mw in molecular_mass_dict.items():
        aa_weight = AASeq.count(amino_acid)*mw
        mass_list.append(aa_weight)
    total_mass = sum(mass_list)
    total_mass = round(total_mass, 2)
    
    return total_mass


def get_biopy_feat(sequence):
    """
    Taken-in dataframe with a column specified as 'SEQUENCE' and uses Biopython to 
    generate additional features from each protein sequence to help better train a machine
    learning model to predict log2FC
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    ProtA = ProteinAnalysis(sequence)
    molwt_biopy = ProtA.molecular_weight()
    aromaticity = ProtA.aromaticity()
    isoelectric_pt = ProtA.isoelectric_point()

    return molwt_biopy, aromaticity, isoelectric_pt


def count_aa_types(sequence):
    '''here, sequence is an element in a list, not a list itself
    Function takes the sequence, counts how many times amino acids occur per each of 4 
    groups and returns values as percentages of the total amino acides in the sequence'''
    l = len(sequence)
    
    ninja = {
        'nonpolar': ['G', 'g', 'A', 'a', 'V', 'v', 'C', 'c', 'P', 'p', 'L', 'l', 'I', 'i', 'M', 'm', 'W', 'w', 'F', 'f'],
        'positive': ['K', 'k', 'R', 'r', 'H', 'h'],
        'polar': ['S', 's', 'T', 't', 'Y', 'y', 'N', 'n', 'Q', 'q'],
        'negative': ['D', 'd', 'E', 'e']
    }
    
    # Initiate lists for counts
    nonpolar_l = []
    positive_l = []
    polar_l = []
    negative_l = []
    
    for i in range(l):
        aa = sequence[i]
        if aa in ninja['nonpolar']:
            nonpolar_l.append(1)
        elif aa in ninja['positive']:
            positive_l.append(1)
        elif aa in ninja['polar']:
            polar_l.append(1)
        else:
            negative_l.append(1)
            
    # Sum the lists to return
    nonpolar_c = sum(nonpolar_l)
    positive_c = sum(positive_l)
    polar_c = sum(polar_l)
    negative_c = sum(negative_l)
    
    # Percentage calculation to return
    nonpolar = nonpolar_c / l *100
    nonpolar = round(nonpolar, 2)
    positive = positive_c / l *100
    positive = round(positive, 2)
    polar = polar_c / l *100
    polar = round(polar, 2)
    negative = negative_c / l *100
    negative = round(negative, 2)
    
    return nonpolar, positive, polar, negative



def get_feat(user_sequence):
    """
    This function takes the sequence that the user is querying and returns a dataframe appended 
    with all of the features used in the predictive model
    """
    import pandas as pd
    mass = get_prot_mass(user_sequence)
    
    molwt_biopy, aromaticity, isoelectric_pt = get_biopy_feat(user_sequence)
    
    nonpolar, positive, polar, negative = count_aa_types(user_sequence)
    
    # Make the features into a Pandas DataFrame
    column_names = ['AA_NP', 'AA_POS', 'AA_POL', 'AA_NEG', 'MW', 'AROM', 'ISO_E']   
    feat_list = [nonpolar, positive, polar, negative, mass, aromaticity, isoelectric_pt]
    feat_df = pd.DataFrame([feat_list], columns=column_names)
    
    return feat_df




def scale_input_feat(X):
    '''
    The function takes in X (our input features), and rescale based on min-max normalization
    it returns the normalized X
    '''
    import pandas as pd
    from sklearn import preprocessing
    #returns a numpy array for X (needed to use the min_max_scaler)
    X_arr = X.values 

    X_col_names = list(X.columns.values.tolist()) #get column names to then put back into X_norm

    #min-max normalization (rescaling) of input features
    min_max_scaler = preprocessing.MinMaxScaler()
    X_scaled = min_max_scaler.fit_transform(X_arr)
    X_norm= pd.DataFrame(X_scaled)

    #put back the original column names
    X_norm.columns = X_col_names
    
    return X_norm


def bagging_regr(X_user):
    '''
    the function takes in the master csv and user input
    csv_file should be inputted as a string = 'compiled_features_complete.csv'
    test ratio, random state and n_estimator are set (from previous ML optimization)
    returning the predicted output based on user input
    '''
    import pandas as pd
    from sklearn import preprocessing
    
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_squared_error, r2_score
    from sklearn.ensemble import BaggingRegressor
    from sklearn.tree import DecisionTreeRegressor
    # Open and load dataset
    bacterial_csv = pd.read_csv('UWInSpace_ModelData.csv')
    df = pd.DataFrame(data=bacterial_csv)
    
    #assign input (X) /output (y) features
    X= df[['AA_NP','AA_POS', 'AA_POL','AA_NEG', 'MW', 'AROM', 'ISO_E']]
    y= df['LOG2FC']
    
    #Scale input features
    X_arr = X.values #returns a numpy array for X (needed to use the min_max_scaler)

    X_col_names = list(X.columns.values.tolist()) #get column names to then put back into X_norm

    #min-max normalization (rescaling) of input features
    min_max_scaler = preprocessing.MinMaxScaler()
    X_scaled = min_max_scaler.fit_transform(X_arr)
    X_norm= pd.DataFrame(X_scaled)

    #put back the original column names
    X_norm.columns = X_col_names
    
  
    #set Bagging regressor parameters, from ML training: 
    test_ratio = 0.30
    seed_random = 42
    n_estim= 20
    
    X_train, X_test, y_train, y_test = train_test_split(X_norm, y, test_size=test_ratio, random_state=seed_random, shuffle=True)
    #Model is Bagging Regressor, base estimator is Decision Tree regressor
    model = BaggingRegressor(base_estimator=DecisionTreeRegressor(),n_estimators=n_estim, random_state=seed_random)
    model.fit(X_train, y_train)
    #y_testpredict = model.predict(X_test)
       
      
    user_predict = model.predict(X_user)
    
    return user_predict


def predict_log2fc(user_input, email):
    """
    This function takes the user input of either a protein accesion number or 
    an amino acid sequence and predicts how much the sequence will change after 
    being in space (near-zero gravity and increased radiation exposure)
    """
    import Bio
    from Bio import Entrez
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import math
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt   

    # Import Scikit-Learn library for the regression model
    import sklearn   
    from sklearn import preprocessing #sklearn for normalization function
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_squared_error, r2_score
    #for Bagging regressor
    from sklearn.ensemble import BaggingRegressor
    from sklearn.tree import DecisionTreeRegressor
    
    # Get user sequence
    user_sequence = usr_seq(user_input, email)
    
    # Get features
    user_features_df = get_feat(user_sequence)
    
    # Predict log2fc
    prediction = bagging_regr(user_features_df)
    prediction = prediction[0]
    prediction = round(prediction, 4)
    
    user_features_df['LOG2FC'] = [prediction]
    
    print('The predicted log2fc is ', prediction)
    
    return prediction, user_features_df



def multi_pred_log2fc(user_input, email):
    """Input the proteins of interest as a list.  This function will iterate through 
    and return a list of predicted log2fc's.
    """
    
    input_length = len(user_input)
    
    # Initiate results list
    result_log2fc = []
    
    for i in range(input_length):
        result, df = predict_log2fc(user_input[i], email)
        result_log2fc.append(result)
    
    return result_log2fc


def cluster_visuals(user_df, user_genecluster):
    import pandas as pd
    import seaborn as sns
    from matplotlib import pyplot as plt

    # set seaborn as default plotting language
    sns.set()

    # function definition
    # NASA_data is the "compiled features complete CSV", maybe need to change to final version of this

    UD = user_df
    NASA_data = pd.read_csv('UWInSpace_ModelData.csv')
    
    ND_ABC = NASA_data[NASA_data['GENENAME'] == 'ABC_transporter']
    ND_CYT = NASA_data[NASA_data['GENENAME'] == 'cytochrome']
    ND_DNA = NASA_data[NASA_data['GENENAME'] == 'DNA_polymerase']
    ND_EFF = NASA_data[NASA_data['GENENAME'] == 'efflux_transporter']
    ND_ELE = NASA_data[NASA_data['GENENAME'] == 'electron_transfer']
    ND_FLA = NASA_data[NASA_data['GENENAME'] == 'flagellar_motility']
    ND_NAD = NASA_data[NASA_data['GENENAME'] == 'NADH_NADPH']
    ND_RNA = NASA_data[NASA_data['GENENAME'] == 'RNA_polymerase']
    
    fig, axes = plt.subplots(1, 7, figsize= (15, 5), sharey=True)
    fig.suptitle('Comparison of features per user-specified gene cluster')
    
    if user_genecluster == 'ABC_transporter':
        ND = ND_ABC
    elif user_genecluster == 'cytochrome':
        ND = ND_CYT
    elif user_genecluster == 'DNA_polymerase':
        ND = ND_DNA
    elif user_genecluster == 'efflux_transporter':
        ND = ND_EFF
    elif user_genecluster == 'electron_transfer':
        ND = ND_ELE
    elif user_genecluster == 'flagellar_motility':
        ND = ND_FLA
    elif user_genecluster == 'NADH_NADPH':
        ND = ND_NAD
    elif user_genecluster == 'RNA_polymerase':
        ND = ND_RNA
    else: 
        print('Please check your user-defined gene cluster and ensure that it matches the possible 8 categories')
    
    # user protein df
    p0a = sns.scatterplot(ax=axes[0], data=UD, x='AA_NP', y='LOG2FC', hue='LOG2FC', palette=['red'], legend=False)
    p1a = sns.scatterplot(ax=axes[1], data=UD, x='AA_POL', y='LOG2FC', hue='LOG2FC', palette=['red'], legend=False)
    p2a = sns.scatterplot(ax=axes[2], data=UD, x='AA_POS', y='LOG2FC', hue='LOG2FC', palette=['red'], legend=False)
    p3a = sns.scatterplot(ax=axes[3], data=UD, x='AA_NEG', y='LOG2FC', hue='LOG2FC', palette=['red'], legend=False)
    p4a = sns.scatterplot(ax=axes[4], data=UD, x='MW', y='LOG2FC', hue='LOG2FC', palette=['red'], legend=False)
    p5a = sns.scatterplot(ax=axes[5], data=UD, x='AROM', y='LOG2FC', hue='LOG2FC', palette=['red'], legend=False)
    p6a = sns.scatterplot(ax=axes[6], data=UD, x='ISO_E', y='LOG2FC', hue='LOG2FC', palette=['red'])
    
    p0 = sns.scatterplot(ax=axes[0], data=ND, x='AA_NP', y='LOG2FC', hue='ORG', legend=False)
    p1 = sns.scatterplot(ax=axes[1], data=ND, x='AA_POL', y='LOG2FC', hue='ORG', legend=False)
    p2 = sns.scatterplot(ax=axes[2], data=ND, x='AA_POS', y='LOG2FC', hue='ORG', legend=False)
    p3 = sns.scatterplot(ax=axes[3], data=ND, x='AA_NEG', y='LOG2FC', hue='ORG', legend=False)
    p4 = sns.scatterplot(ax=axes[4], data=ND, x='MW', y='LOG2FC', hue='ORG', legend=False)
    p5 = sns.scatterplot(ax=axes[5], data=ND, x='AROM', y='LOG2FC', hue='ORG', legend=False)
    p6 = sns.scatterplot(ax=axes[6], data=ND, x='ISO_E', y='LOG2FC', hue='ORG')
    
    return
#################################################################################

def test_get_seq():
    seq = get_seq("karanjia@uw.edu", 728886557)
    expect = 'mlrsmltasttlnqlqqqidtissnlsnsnttgykakdtnfselvrqqfdqvdekneevakarktppglrlgvgammssrlvsdqgsiqktdrdldiaftspyqylqvnvngnrqytrdgalyvtpsaananqlqlvtgngypvldengntvnidssmknitinkngtltasdgnavqrfnlgvvqvnnpqelksegnnlfsidnaaafeelnganrqnigmqqgslemsnvdiseqmtdlitsqrsyqlnsrtitmgdqmlglinsvr'
    assert seq == expect, 'the get_seq function is not working'
    return


def test_usr_seq():
    # Case 1: input is EID
    email = 'jaking11@uw.edu'
    usr_in = 15599626
    result = usr_seq(usr_in, email)
    assert result == 'mnkfmawvdarfpatkmwedhlskyyapknfnfwyffgslallvlvnqiltgiwltmsftpsaeeafasveyimrdvdygwiirymhstgasaffivvylhmfrgllygsyqkprelvwifgmliylalmaeafmgyllpwgqmsywgaqviislfgaipvvgedlaqwirgdflisgitlnrffalhvialpivllglvvlhilalhevgsnnpdgvdikkkkdengvpldgiafhpyytvkdivgvvvflfifctvifffpemggyflekpnfemanqfktpehiapvwyftpfyailravpdklmgvvamgaaiavlfvlpwldrspvrsirykgwlsklwlvifavsfvilgyygaqapsplgttlsrvctvlyfaffilmpfytrmektkpvpervtg', 'the usr_seq function does not work'
    return


def test_get_prot_mass():
    import math
    seq = 'TACO'
    result = get_prot_mass(seq)
    expect = 512.62
    assert math.isclose(result, expect), 'the get_prot_mass function does not work'
    return


def test_get_biopy_feat():
    import math
    seq = 'TACO'
    molwt_biopy, aromaticity, isoelectric_pt = get_biopy_feat(seq)
    expect = 530.6381
    assert math.isclose(molwt_biopy, expect), 'the get_biopy_feat function does not work'
    return


def test_count_aa_types():
    import math
    seq = 'TACO'
    nonpolar, positive, polar, negative = count_aa_types(seq)
    expect = 25.0
    assert math.isclose(polar, expect), 'the count_aa_types function does not work'
    return


def test_get_feat():
    import math
    seq = 'TACO'
    df = get_feat(seq)
    result = df['MW']
    expect = 512.62
    assert math.isclose(result, expect), 'the get_feat function is not working'
    return


def test_scale_input_feat():
    import math
    import numpy as np
    import pandas as pd

    arr = np.array([[1, 2], [2, 3], [3, 1]])
    df = pd.DataFrame(arr)
    df_norm = scale_input_feat(df)    
    result = df_norm[0].max()
    expect = 1.0
    assert math.isclose(result, expect), 'the scale_input_feat function is not working'
    return


def test_predict_log2fc():
    import math
    user_input = 15600267
    email = 'jaking11@uw.edu'
    prediction, user_features_df = predict_log2fc(user_input, email)
    expect = 0.168
    assert math.isclose(prediction, expect), 'the predict_log2fc function does not work'
    return
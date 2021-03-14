# This is a file of test functions for all of the functions in our class project.

#################################################################################
# NEED TO COPY IN utils.py
#################################################################################

def test_usr_seq():
    # Case 1: input is EID
    usr_in = 15599626
    result = usr_seq(usr_in)
    assert result is str, 'the usr_seq function does not work'
    # Case 2: input is sequence
    usr_in = 'juliaisthecoolest'
    result = usr_seq(usr_in)
    assert result == usr_in, 'the usr_seq function does not work'
    return
    
def test_get_biopy_feat():
    seq = 'mlpamrtgllcallgvtapawaeyvtvisfggankeaqetafykpfksatgnrvvhgsyngdlaklkrmveishvswdvveveapelargceeglfekldmakvgdpadfvpgavqpcgvgifvwttllaynpgkvagspqgwadfwdvkkfpgkrglrwgakyslefalmadgvapkdvyqtlatpagverafrkldelkpyihwwksgqdpvrdladgtvvmssayngriaaaqaekqrlamvwsggvydfdfwalpvgvwkkqlaeefirfasqpeqqkafaeniaygpanrkavglldpqvaanlptapqnmqnavgmnvafwaehgealeqrfqnwakr'
    molwt_biopy, aromaticity, isoelctric_pt = get_biopy_feat(seq)
    #print(molwt_biopy)
    #print(aromaticity)
    #print(isoelctric_pt)
    expect_molwt = 37832.81
    assert math.isclose(molwt_biopy, expect_molwt, rel_tol=0.05), 'the get_biopy_feat function does not work'
    return

def test_aa_to_df():
    # Initiate dataframe for test function
    fasta = ['ALKNRRTV', 'AQWWWQP', 'DFITNML', 'CIYRREMWLNLATD', 'FDTWMNQI', 
             'PMPSKSL', 'TFGSRECVPDSLIE', 'CATS', 'GEWLIA', 'A', 'KING']
    df = pd.DataFrame(fasta, columns=['SEQUENCE'])
    
    # Test the function
    df = aa_to_df(df)
    result = df['AA_NEG'][0]
    expect = 0
    assert math.isclose(result, expect), 'The test_aa_to_df function does not work properly.'
    return

def test_scale_input_feat():
    arr = [[3, 4, 2], [1, 7, 0], [6, 5, 2]]
    names = ['1st', '2nd', '3rd']
    df = pd.DataFrame(arr, columns=names)
    scaled_df = scale_input_feat(df)
    result = scaled_df['1st'][0]
    expect = 0.4
    assert math.isclose(result, expect), 'the scale_input_feat function does not work'
    return
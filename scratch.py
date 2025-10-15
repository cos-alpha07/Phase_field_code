import pickle, os
from pandas import DataFrame
from MAIN import *


def save_my_results(U1, U2, Phi, his_par, srmat, tsmat, stp):
    
    gpts = get_gpts_crds() #recheck
    srmat_df = DataFrame(srmat, columns=['disp', 'load'])
    tsmat_df = DataFrame(tsmat, columns=['step', 'iterations', 'step time'])
    
                                                             
    u1[f"Step-{step_num}"] = pd.Series(U[[dim * i for i in range(node_count)]],
                                       index=range(len(node_mat)))
    if dim > 1:
        u2[f"Step-{step_num}"] = pd.Series(U[[dim * i + 1 for i in range(node_count]],
                                           index=range(len(node_mat)))
        
    dmg[f"Step-{step_num}"] = pd.Series(Phi, index=range(node_count))
    
    cys = ["X", "Y", "Z"]
    nds = [f"n{i+1}" for i in range(el_data.shape[1] - 1)]
    node_df = pd.DataFrame(nd_data[:, : dim + 1], columns=["Node"] + cys[:dim])
    elem_df = pd.DataFrame(el_data, columns=list(["Element"] + nds))
    gpt_df = pd.DataFrame(gpts, columns=[f"{cys[i]}coord" for i in range(dim)], index=range(len(gpts)))
    
    # dumping data-pickles
    dumps = dict(
        [
            ("u1", locals()["u1"]),
            ("u2", locals()["u2"]),
            ("SRMat_df", locals()["srmat_df"]),
            ("TSMat_df", locals()["tsmat_df"]),
            ("node_df", locals()["node_df"]),
            ("elem_df", locals()["elem_df"]),
            ("gpt_df", locals()["gpt_df"]),
            ("dmg", locals()["dmg"]),
            ("stepHistory", locals()["stepHistory"]),
        ]
    )
    # Define the file path where the selected variables will be saved
    dump_file_path = rf"{H:\kp\Phase_field_code\results}\{probname}\data.pkl"
    
    # Save the selected variables to a file
    with open(dump_file_path, "wb") as f:
        pickle.dump(dumps, f)
        
    print("data saved\n")


    

    
    
    





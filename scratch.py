import pickle
import os


def save_my_results(U, Phi, his_par, srmat, tsmat, stp):
    folder = 'saved_data'
    os.makedirs(folder, exist_ok=True)
    
    data_obj = {f'displacement data at step-{stp}': U,
            f'PF variable data at step-{stp}': Phi,
            f'History parameter variable at step-{stp}': his_par,
            f'SR data at step-{stp}': srmat,
            f'TS data at step-{stp}': tsmat}
    
    for filename, data in data_obj.items():
        file_path = os.path.join(folder, f'{filename}.pkl')
        pickle.dump(data, open(file_path, 'wb'))
    

    
    
    





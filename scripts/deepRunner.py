import sys
import argparse
from DeepTCR.DeepTCR import DeepTCR_U
from DeepTCR.DeepTCR import DeepTCR_SS
from DeepTCR.DeepTCR import DeepTCR_WF

def DeepTCR_supervised(dpath, study_id):
    # Instantiate training object
    deeptcr = DeepTCR_WF(study_id)

    # Load TCR Data from directories
    deeptcr.Get_Data(directory=dpath, Load_Prev_Data=False,
        aa_column_beta=1,count_column=2,v_beta_column=7,d_beta_column=14,
        j_beta_column=21,type_of_data_cut='Num_Seq',
        data_cut = 5000)
    # Train and perform 10-Fold validation
    folds = 10
    size_of_net = 'small'
    num_concepts=64
    hinge_loss_t = 0.1
    train_loss_min=0.1
    deeptcr.K_Fold_CrossVal(combine_train_valid=True, 
        hinge_loss_t = hinge_loss_t,
        train_loss_min = train_loss_min,
        folds=folds,
        num_concepts=num_concepts, size_of_net=size_of_net)

    deeptcr.Representative_Sequences(top_seq=100, motif_seq = 100)
    deeptcr.HeatMap_Sequences()
    deeptcr.HeatMap_Samples()
    deeptcr.Repertoire_Dendrogram()
    deeptcr.UMAP_Plot()
    deeptcr.UMAP_Plot_Samples()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='DeepTCR Runner')
    parser.add_argument('-d', '--dpath', type = str, help='Path to DeepTCR formatted data')
    parser.add_argument('-s', '--study_id', type = str, help='Study ID. Will be the name of the output folder created by DeepTCR')
    parser.add_argument('-m', '--mode', type = str, help='foo help', default='supervised')
    args = parser.parse_args()

    if args.mode == "supervised":
        DeepTCR_supervised(args.dpath, args.study_id)
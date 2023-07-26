#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 15:34:42 2019

@author: yanying
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import argparse
import itertools
import os
import time 
import seaborn as sns
import pandas
import sklearn.model_selection
import sklearn.metrics
from sklearn import linear_model
from scipy.stats import spearmanr,pearsonr
from collections import defaultdict
import shap
import sys
from sklearn.preprocessing import LabelEncoder,StandardScaler
import pickle
from sklearn.feature_selection import VarianceThreshold
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 400
mpl.rcParams['font.sans-serif']='Arial'
mpl.rcParams['font.size']=14
mpl.rcParams['legend.title_fontsize']=10
mpl.rcParams['legend.fontsize']=10
mpl.rcParams['xtick.labelsize']=12
mpl.rcParams['ytick.labelsize']=12
import warnings
warnings.filterwarnings('ignore')
start_time=time.time()
nts=['A','T','C','G']
items=list(itertools.product(nts,repeat=2))
dinucleotides=list(map(lambda x: x[0]+x[1],items))
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
        
parser = MyParser(usage='python %(prog)s datasets [options]',formatter_class=argparse.RawTextHelpFormatter,description="""
This is used to train optimized models from auto-sklearn and other model types with individual or fused datasets, and evaluate with 10-fold cross-validation. 

Example: python BE_kill_machine_learning.py -o test -r regressor.pkl
                  """)
parser.add_argument("data_csv", help="data csv file")
parser.add_argument("-o", "--output", default="results", help="output folder name. default: results")
parser.add_argument("-f","--folds", type=int, default=10, help="Fold of cross validation, default: 10")
parser.add_argument("-t","--test_size", type=float, default=0.2, help="Test size for spliting datasets, default: 0.2")
parser.add_argument("-s","--seq_choice", type=str, default=None, help="If adding killing guide sequence features, if yes, use -s add_kill. default: None")
parser.add_argument("-r","--regressor", type=str, default=None, help="Saved regressor from autosklearn, default: None")


args = parser.parse_args()
output_file_name = args.output
folds=args.folds
test_size=args.test_size
regressor=args.regressor
data_csv=args.data_csv
seq_choice=args.seq_choice

try:
    os.mkdir(output_file_name)
except:
    overwrite=input("File exists, do you want to overwrite? (y/n)")
    if overwrite == "y":
        os.system("rm -r %s"%output_file_name)
        os.mkdir(output_file_name)
    elif overwrite =="n":
        output_file_name=input("Please give a new output file name:")
        os.mkdir(output_file_name)
    else:
        print("Please input valid choice..\nAbort.")
        sys.exit()
def self_encode(sequence):#one-hot encoding for single nucleotide features
    integer_encoded=np.zeros([len(sequence),4],dtype=np.float64)
    nts=['A','T','C','G']
    for i in range(len(sequence)):
        integer_encoded[i,nts.index(sequence[i])]=1
    sequence_one_hot_encoded = integer_encoded.flatten()
    return sequence_one_hot_encoded

def dinucleotide(sequence):#encoding for dinucleotide features
    nts=['A','T','C','G']
    items=list(itertools.product(nts,repeat=2))
    dinucleotides=list(map(lambda x: x[0]+x[1],items))
    encoded=np.zeros([(len(nts)**2)*(len(sequence)-1)],dtype=np.float64)
    for nt in range(len(sequence)-1):
        if sequence[nt] == 'N' or sequence[nt+1] =='N':
            print(sequence)
            continue
        encoded[nt*len(nts)**2+dinucleotides.index(sequence[nt]+sequence[nt+1])]=1
    return encoded

def DataFrame_input(df,coding_strand=1):
    ###keep guides for essential genes
    logging_file= open(output_file_name + '/log.txt','a')
    print(df.shape)
    df=df[(df['essentiality_M9']==1)|(df['essentiality_MOPS']==1)]
    print(df.shape)
    sequences=list(dict.fromkeys(df['sgRNA']))
     ### one hot encoded sequence features
    PAM_encoded=[]
    sequence_encoded=[]
    dinucleotide_encoded=[]
    for i in df.index:
        PAM_encoded.append(self_encode(df['PAM'][i]))
        sequence_encoded.append(self_encode(df['sgRNA'][i]))
        dinucleotide_encoded.append(dinucleotide(df['sequence_30nt'][i]))
        
        if df['essentiality_M9'][i]==1 and df['essentiality_MOPS'][i]==0:
            df.at[i,'logFC']=df['M9-killing.logFC'][i]
        elif df['essentiality_M9'][i]==0 and df['essentiality_MOPS'][i]==1:
            df.at[i,'logFC']=df['MOPS-killing.logFC'][i]
        elif df['essentiality_M9'][i]==1 and df['essentiality_MOPS'][i]==1:
            # print('error')
            df.at[i,'logFC']=np.mean([df['MOPS-killing.logFC'][i],df['M9-killing.logFC'][i]])
        
        df.at[i,'geneid']=int(df['geneid'][i][1:])
        df.at[i,'guideid']=sequences.index(df['sgRNA'][i])
        
    #check if the length of gRNA and PAM from all samples is the same
    if len(list(set(map(len,list(df['PAM'])))))==1:
        PAM_len=int(list(set(map(len,list(df['PAM']))))[0])
    else:
        print("error: PAM len")
    if len(list(set(map(len,list(df['sgRNA'])))))==1:   
        sequence_len=int(list(set(map(len,list(df['sgRNA']))))[0])
    else:
        print("error: sequence len")
    if len(list(set(map(len,list(df['sequence_30nt'])))))==1:   
        dinucleotide_len=int(list(set(map(len,list(df['sequence_30nt']))))[0])
    else:
        print("error: sequence len")
    guideids=np.array(list(df['guideid']))
    
    df=df.dropna(subset=['logFC'])
    print(df.shape)
    logging_file.write("Number of selected guides: %s\n" % df.shape[0])
    
    y=np.array(df['logFC'],dtype=float)
    
    headers=list()
    ### add one-hot encoded sequence features columns
    PAM_encoded=np.array(PAM_encoded)
    sequence_encoded=np.array(sequence_encoded)
    dinucleotide_encoded=np.array(dinucleotide_encoded)
    X=np.c_[sequence_encoded,PAM_encoded,dinucleotide_encoded]
    ###add one-hot encoded sequence features to headers
    for i in range(sequence_len):
        for j in range(len(nts)):
            headers.append('sequence_%s_%s'%(i+1,nts[j]))
    for i in range(PAM_len):
        for j in range(len(nts)):
            headers.append('PAM_%s_%s'%(i+1,nts[j]))
    for i in range(dinucleotide_len-1):
        for dint in dinucleotides:
            headers.append(dint+str(i+1)+str(i+2))
    X=pandas.DataFrame(data=X,columns=headers,dtype=float)
    logging_file.write("Number of guides for training: %s \n" % X.shape[0])
    logging_file.write("Number of features: %s\n" % len(headers))
    logging_file.write("Features: "+",".join(headers)+"\n\n")
    return X, y, headers, guideids


def Evaluation(output_file_name,y,predictions,name):
    #scores
    spearman_rho,spearman_p_value=spearmanr(y, predictions)
    pearson_rho,pearson_p_value=pearsonr(y, predictions)
    y=np.array(y)
    
    # scatter plot
    plt.figure() 
    sns.set_palette("PuBu",2)
    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_xDist = plt.subplot(gs[0, :2],sharex=ax_main)
    ax_yDist = plt.subplot(gs[1:3, 2],sharey=ax_main)
    ax_main.scatter(y,predictions,edgecolors='white',alpha=0.8)
    ax_main.set(xlabel='Experimental log2FC',ylabel='Predicted log2FC')
    ax_xDist.hist(y,bins=70,align='mid',alpha=0.7)
    ax_xDist.set(ylabel='count')
    ax_xDist.tick_params(labelsize=6,pad=2)
    ax_yDist.hist(predictions,bins=70,orientation='horizontal',align='mid',alpha=0.7)
    ax_yDist.set(xlabel='count')
    ax_yDist.tick_params(labelsize=6,pad=2)
    ax_main.text(0.55,0.03,"Spearman R: {0}".format(round(spearman_rho,2)),transform=ax_main.transAxes,fontsize=10)
    ax_main.text(0.55,0.10,"Pearson R: {0}".format(round(pearson_rho,2)),transform=ax_main.transAxes,fontsize=10)
    plt.savefig(output_file_name+'/'+name+'_scatterplot.png',dpi=300)
    plt.close()
    

def SHAP(estimator,X,headers):
    X=pandas.DataFrame(X,columns=headers)
    X.to_csv(output_file_name+"/shap_samples.csv",sep='\t',index=False)
    X=X.astype(float)
    explainer=shap.TreeExplainer(estimator)
    shap_values = explainer.shap_values(X,check_additivity=False)
    values=pandas.DataFrame({'shap_values':np.mean(np.absolute(shap_values),axis=0),'features':headers})
    values.to_csv(output_file_name+"/shap_value_mean.csv",index=False,sep='\t')
    
    shap.summary_plot(shap_values, X, plot_type="bar",show=False,color_bar=True,max_display=10)
    plt.subplots_adjust(left=0.35, top=0.95)
    plt.savefig(output_file_name+"/shap_value_bar.svg",dpi=400)
    plt.close()
    
    for i in [10,15,30]:
        shap.summary_plot(shap_values, X,show=False,max_display=i,alpha=0.5)
        plt.subplots_adjust(left=0.45, top=0.95,bottom=0.2)
        plt.yticks(fontsize='medium')
        plt.xticks(fontsize='medium')
        plt.savefig(output_file_name+"/shap_value_top%s.svg"%(i),dpi=400)
        plt.savefig(output_file_name+"/shap_value_top%s.png"%(i),dpi=400)
        plt.close()    
    
    shap_values=pandas.DataFrame(shap_values,columns=headers)
    shap_values.to_csv(output_file_name+"/shap_values.csv",sep='\t',index=False)
    

def main():
    open(output_file_name + '/log.txt','a').write("Python script: %s\n"%sys.argv[0])
    open(output_file_name + '/log.txt','a').write("Parsed arguments: %s\n\n"%args)
    df1=pandas.read_csv(data_csv,sep="\t")
    open(output_file_name + '/log.txt','a').write("Total number of guides in dataset %s: %s\n"% (data_csv,df1.shape[0]))
    training_df = df1.sample(frac=1,random_state=np.random.seed(111)).reset_index(drop=True)
    #dropping unnecessary features and encode sequence features
    X,y,headers,guideids=DataFrame_input(training_df)
    open(output_file_name + '/log.txt','a').write("Data input Time: %s seconds\n\n" %('{:.2f}'.format(time.time()-start_time)))  
    
    numerical_indicator=["gene_GC_content","distance_operon","distance_operon_perc","operon_downstream_genes","ess_gene_operon","gene_length","gene_expression_min","gene_expression_max",\
                          'distance_start_codon','distance_start_codon_perc','guide_GC_content','MFE_hybrid_seed','MFE_homodimer_guide','MFE_hybrid_full','MFE_monomer_guide',\
                        'homopolymers']
    dtypes=dict()
    for feature in headers:
        if feature not in numerical_indicator:
            dtypes.update({feature:int})
    X=pandas.DataFrame(data=X,columns=headers)
    X=X.astype(dtypes)
    
    ##optimized models from auto-sklearn
    if regressor ==None:
        from sklearn.experimental import enable_hist_gradient_boosting
        from sklearn.ensemble import HistGradientBoostingRegressor
        from sklearn.ensemble import RandomForestRegressor
        estimator=pickle.load(open(regressor,'rb'))
        print(estimator.get_params())
        params=estimator.get_params()
        
        params.update({"random_state":np.random.seed(111)})
        if 'max_iter' in params.keys():
            params.update({'max_iter':512})
        if 'early_stop' in params.keys():
            params.pop('early_stop', None)
        if 'max_depth' in params.keys():
            if params['max_depth']=='None':
                params['max_depth']=None
        if 'max_leaf_nodes' in params.keys():
            if params['max_leaf_nodes']=='None':
                params['max_leaf_nodes']=None
        if 'Gradient Boosting' in str(estimator):
            from sklearn.experimental import enable_hist_gradient_boosting
            from sklearn.ensemble import HistGradientBoostingRegressor
            estimator=HistGradientBoostingRegressor(**params)
        elif 'Extra Trees' in str(estimator):
            estimator=sklearn.ensemble.ExtraTreesRegressor(**params)
        elif 'Random Forest' in str(estimator):
            from sklearn.ensemble import RandomForestRegressor
            estimator=RandomForestRegressor(**params)
            
        open(output_file_name + '/log.txt','a').write("regressor:"+str(estimator)+"\n")
    scaler=StandardScaler()
    selector = VarianceThreshold()
    open(output_file_name + '/log.txt','a').write("selector:"+str(selector)+"\n")
    open(output_file_name + '/log.txt','a').write("scaler:"+str(scaler)+"\n")
    
        
    open(output_file_name + '/log.txt','a').write("Estimator:"+str(estimator)+"\n")
    X_df=pandas.DataFrame(data=np.c_[X,y,guideids],columns=headers+['log2FC','guideid'])
    #k-fold cross validation
    evaluations=defaultdict(list)
    kf=sklearn.model_selection.KFold(n_splits=folds, shuffle=True, random_state=np.random.seed(111))
    guideid_set=list(set(guideids))
    for train_index, test_index in kf.split(guideid_set):##split the combined training set into train and test based on guideid
        train_index=np.array(guideid_set)[train_index]
        test_index=np.array(guideid_set)[test_index]
        X_train = X_df[X_df['guideid'].isin(train_index)]
        y_train=X_train['log2FC']
        X_train=X_train[headers]
        X_train=X_train.astype(dtypes)
        
        X_train=selector.fit_transform(X_train)
        X_train=scaler.fit_transform(X_train)
        
        test = X_df[X_df['guideid'].isin(test_index)]
        X_test=test.copy()
        y_test=X_test['log2FC']
        X_test=X_test[headers]
        X_test=X_test.astype(dtypes)
        X_test=selector.transform(X_test)
        X_test=scaler.transform(X_test)
    
        estimator = estimator.fit(np.array(X_train,dtype=float),np.array(y_train,dtype=float))
        predictions = estimator.predict(np.array(X_test,dtype=float))
        # print(spearmanr(y_test, predictions,nan_policy='omit')[0])
        evaluations['Rs'].append(spearmanr(y_test, predictions)[0])
        # Evaluation(output_file_name,y_test,predictions,"X_test_kfold")
     
    evaluations=pandas.DataFrame.from_dict(evaluations)
    evaluations.to_csv(output_file_name+'/iteration_scores.csv',sep='\t',index=True)
    
    logging_file= open(output_file_name + '/log.txt','a')
    ##split the combined training set into train and test
    ### save model trained with all guides
    X_all=X_df[headers]
    y_all=X_df['log2FC']
    # X_all.to_csv(output_file_name+"/shap_samples_nonscaled.csv",sep='\t',index=False)
    X_all=selector.fit_transform(X_all)
    mask=selector.get_support()
    if False not in mask:
        new_headers=headers
    else:
        if len(mask)==len(headers):
            new_headers=[]
            for i in range(len(mask)):
                if mask[i]:
                    new_headers.append(headers[i])
    print(len(new_headers))
    X_all=scaler.fit_transform(X_all)
    estimator = estimator.fit(np.array(X_all,dtype=float),np.array(y_all,dtype=float))
    if os.path.isdir(output_file_name+'/saved_model')==False:  
        os.mkdir(output_file_name+'/saved_model')
    pickle.dump(selector, open(output_file_name+'/saved_model/selector.sav', 'wb'))
    pickle.dump(scaler, open(output_file_name+'/saved_model/scaler.sav', 'wb'))
    pickle.dump(estimator, open(output_file_name+'/saved_model/estimator.sav', 'wb'))
    pickle.dump(headers, open(output_file_name+'/saved_model/headers.sav', 'wb'))
    
    SHAP(estimator,X_all,new_headers)
    Evaluation(output_file_name,y_all,estimator.predict(np.array(X_all,dtype=float)),"X_all")
    logging_file.close()
    

if __name__ == '__main__':
    main()
    open(output_file_name + '/log.txt','a').write("Execution Time: %s seconds\n" %('{:.2f}'.format(time.time()-start_time)))    
#%%

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
import logging
import pandas
import sys
import sklearn.model_selection
import sklearn.metrics
import autosklearn.regression
import autosklearn.classification
import autosklearn.metrics
from sklearn.preprocessing import LabelEncoder,StandardScaler
from sklearn.feature_selection import VarianceThreshold
import pickle
from scipy.stats import spearmanr,pearsonr
start_time=time.time()
import warnings
warnings.filterwarnings('ignore')
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
        
parser = MyParser(usage='python %(prog)s datasets [options]',formatter_class=argparse.RawTextHelpFormatter,description="""
This is used to optimize models using auto-sklearn (tested version 0.14.6).

ensemble_size, folds, per_run_time_limit, time_left_for_this_task, include_estimators, and include_preprocessors are parameters for auto-sklearn. More description please check the API of auto-sklearn (https://automl.github.io/auto-sklearn/master/api.html)

Example: python conditional_essentiality_screen_autosklearn.py gRNAs.csv
                  """)
parser.add_argument("data_csv", help="data csv file")
parser.add_argument("-o", "--output", default="results", help="output folder name. default: results")
parser.add_argument("-e","--ensemble_size", type=int, default=1, help="Ensemble size, default: 1")
parser.add_argument("-f","--folds", type=int, default=10, help="Fold of cross validation, default: 10")
parser.add_argument("-prt","--per_run_time_limit", type=int, default=360, help="per_run_time_limit (in second), default: 360")
parser.add_argument("-ptt","--time_left_for_this_task", type=int, default=3600, help="time_left_for_this_task (in second), default: 3600")
parser.add_argument("-inest","--include_estimators", type=str, default=None, help="estimators to be included in auto-sklearn. Multiple input separated by ','. If None, then include all. Default: None")
parser.add_argument("-inprepro","--include_preprocessors", type=str, default=None, help="preprocessors to be included in auto-sklearn. Multiple input separated by ','. If None, then include all. Default: None")
args = parser.parse_args()
data_csv=args.data_csv
output_file_name = args.output
folds=args.folds
ensemble_size=args.ensemble_size
per_run_time_limit=args.per_run_time_limit
time_left_for_this_task=args.time_left_for_this_task
include_estimators=args.include_estimators
include_preprocessors=args.include_preprocessors
nts=['A','T','C','G']
### esitmator and preprocessor setting for auto sklearn
if include_estimators != None:
    include_estimators=include_estimators.split(',')
if include_preprocessors != None:
    include_preprocessors=include_preprocessors.split(',')

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


def self_encode(sequence):
    integer_encoded=np.zeros([len(sequence),4],dtype=np.float64)
    nts=['A','T','C','G']
    for i in range(len(sequence)):
        integer_encoded[i,nts.index(sequence[i])]=1
    sequence_one_hot_encoded = integer_encoded.flatten()
    return sequence_one_hot_encoded

def dinucleotide(sequence):
    nts=['A','T','C','G']
    items=list(itertools.product(nts,repeat=2))
    dinucleotides=list(map(lambda x: x[0]+x[1],items))
    encoded=np.zeros([(len(nts)**2)*(len(sequence)-1)],dtype=np.float64)
    for nt in range(len(sequence)-1):
        encoded[nt*len(nts)**2+dinucleotides.index(sequence[nt]+sequence[nt+1])]=1
    return encoded


def DataFrame_input(df):
    logging_file= open(output_file_name + '/log.txt','a')
    ###keep essential genes
    print(df.shape)
    df=df[(df['essentiality_M9']==1)|(df['essentiality_MOPS']==1)]
    print(df.shape)
    
    ### adding
    PAM_encoded=[]
    sequence_encoded=[]
    dinucleotide_encoded=[]
    kill_sequence_encoded=[]
    sequences=list(dict.fromkeys(df['sgRNA']))
    nts=['A','T','C','G']
    for i in df.index:
        PAM_encoded.append(self_encode(df['PAM'][i]))
        sequence_encoded.append(self_encode(df['sgRNA'][i]))
        dinucleotide_encoded.append(dinucleotide(df['sequence_30nt'][i]))
        
        if df['essentiality_M9'][i]==1 and df['essentiality_MOPS'][i]==0:
            df.at[i,'logFC']=df['M9-killing.logFC'][i]
        elif df['essentiality_M9'][i]==0 and df['essentiality_MOPS'][i]==1:
            df.at[i,'logFC']=df['MOPS-killing.logFC'][i]
        elif df['essentiality_M9'][i]==1 and df['essentiality_MOPS'][i]==1:
            df.at[i,'logFC']=np.mean([df['MOPS-killing.logFC'][i],df['M9-killing.logFC'][i]])
        df.at[i,'geneid']=int(df['geneid'][i][1:])
        df.at[i,'guideid']=sequences.index(df['sgRNA'][i])
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
    df=df.dropna(subset=['logFC'])
    print(df.shape)
    logging_file.write("Number of selected guides: %s\n" % df.shape[0])
    
    y=np.array(df['logFC'],dtype=float)
    guideids=np.array(list(df['guideid']))
    # 
    ### feat_type for auto sklearn
    feat_type=[]
    ### add sequence columns
    PAM_encoded=np.array(PAM_encoded)
    sequence_encoded=np.array(sequence_encoded)
    dinucleotide_encoded=np.array(dinucleotide_encoded)
    X=np.c_[sequence_encoded,PAM_encoded,dinucleotide_encoded]
    nts=['A','T','C','G']
    headers=list()
    for i in range(sequence_len):
        for j in range(len(nts)):
            headers.append('sequence_%s_%s'%(i+1,nts[j]))
    for i in range(PAM_len):
        for j in range(len(nts)):
            headers.append('PAM_%s_%s'%(i+1,nts[j]))
    items=list(itertools.product(nts,repeat=2))
    dinucleotides=list(map(lambda x: x[0]+x[1],items))
    for i in range(dinucleotide_len-1):
        for dint in dinucleotides:
            headers.append(dint+str(i+1)+str(i+2))
    cat_length=PAM_len*4+sequence_len*4+(dinucleotide_len-1)*4*4
    for i in range(cat_length):
        feat_type.append('Categorical')
        
    X=pandas.DataFrame(X,columns=headers)      
    logging_file.write("Number of features: %s\n" % len(headers))
    logging_file.write("Features: "+",".join(headers)+"\n\n")
    return X, y, guideids,feat_type, headers

def correlation_plot(df,output_file_name):
    plt.figure(figsize=(12,10))
    cor = df.corr()
    sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
    plt.savefig( output_file_name+ '/correlation.png',dpi=400)


def Evaluation(output_file_name,y,predictions,name):
    #scores
    output=open(output_file_name+"/result.txt","a")
    spearman_rho,spearman_p_value=spearmanr(y, predictions)
    pearson_rho,pearson_p_value=pearsonr(y, predictions)
    output.write(name+"\n")
    output.write("spearman correlation rho: "+str(spearman_rho)+"\n")
    output.write("spearman correlation p value: "+str(spearman_p_value)+"\n")
    output.write("pearson correlation rho: "+str(pearson_rho)+"\n")
    output.write("pearson correlation p value: "+str(pearson_p_value)+"\n")
    output.write("r2: "+str(sklearn.metrics.r2_score(y,predictions))+"\n")
    output.write("explained_variance_score score of "+name+" :"+str(sklearn.metrics.explained_variance_score(y, predictions))+"\n")
    output.write("Mean absolute error regression loss score of "+name+" :"+str(sklearn.metrics.mean_absolute_error(y, predictions))+"\n")
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
    


def main():
    open(output_file_name + '/log.txt','a').write("Python script: %s\n"%sys.argv[0])
    open(output_file_name + '/log.txt','a').write("Parsed arguments: %s\n\n"%args)
    training_df=pandas.read_csv(data_csv,sep="\t")
    training_df = training_df.sample(frac=1,random_state=np.random.seed(111)).reset_index(drop=True)
    open(output_file_name + '/log.txt','a').write("training input shape: %s\n"% str(training_df.shape))
    
    X,y,guideids,feat_type,headers=DataFrame_input(training_df)
    open(output_file_name + '/log.txt','a').write("Data input Time: %s seconds\n" %('{:.2f}'.format(time.time()-start_time)))  
    
    X_df=pandas.DataFrame(data=np.c_[X,y,guideids],columns=headers+['logFC','guideid'])
    
    scaler=StandardScaler()
    selector = VarianceThreshold()
    guideid_set=list(set(guideids))
    ##split the combined training set into train and test
    guide_train, guide_test = sklearn.model_selection.train_test_split(guideid_set, test_size=0.2,random_state=np.random.seed(111))  
    X_train = X_df[X_df['guideid'].isin(guide_train)]
    y_train=np.array(X_train['logFC'],dtype=float)
    X_train = X_train[headers]
    X_train=np.array(X_train,dtype=float)
    X_train=selector.fit_transform(X_train)
    X_train=scaler.fit_transform(X_train)
    mask=selector.get_support()
    if False not in mask:
        new_headers=headers
    else:
        if len(mask)==len(headers):
            new_headers=[]
            feat_type=[]
            for i in range(len(mask)):
                if mask[i]:
                    new_headers.append(headers[i])
                    feat_type.append("Categorical")
        print([i for i in headers if i not in new_headers])
        # sys.exit()
        open(output_file_name + '/log.txt','a').write("Number of Features after selector: {0}\n".format(len(new_headers)))
        open(output_file_name + '/log.txt','a').write("Features after selector: "+",".join(new_headers)+"\n\n")
        pickle.dump(new_headers, open(output_file_name+"/headers.pkl", 'wb'))
    
    test = X_df[X_df['guideid'].isin(guide_test)]
    y_test=np.array(test['logFC'],dtype=float)
    X_test = test[headers]
    X_test=np.array(X_test,dtype=float)
    X_test=selector.transform(X_test)
    X_test=scaler.transform(X_test)
    
    # training
    estimator = autosklearn.regression.AutoSklearnRegressor(
            ensemble_size=ensemble_size,
            time_left_for_this_task=time_left_for_this_task,
            per_run_time_limit=per_run_time_limit,
            include = {'feature_preprocessor': ["no_preprocessing"]},
            resampling_strategy='cv',
            resampling_strategy_arguments={'folds': folds},
            tmp_folder=output_file_name+'/autosklearn_regression_example_tmp',
            # output_folder=output_file_name+'/autosklearn_regression_example_out',
            delete_tmp_folder_after_terminate=False,
            # delete_output_folder_after_terminate=False,
            disable_evaluator_output=False,
            memory_limit=None,
            ensemble_nbest=50, seed = 1,
            get_smac_object_callback=None,
            initial_configurations_via_metalearning=25,
            logging_config=None, metadata_directory = None,
            n_jobs= None, smac_scenario_args= None,metric=autosklearn.metrics.r2)
    estimator.fit(X_train.copy(), y_train.copy(),feat_type=feat_type)
    estimator.fit_ensemble(y_train.copy(), task=None, precision=32, dataset_name=None, ensemble_nbest=None, ensemble_size=ensemble_size)
    open(output_file_name + '/log.txt','a').write("Get parameters:\n"+str(estimator.get_params())+"\n\n")
    open(output_file_name + '/log.txt','a').write("Show models: \n"+str(estimator.show_models())+"\n\n")
    open(output_file_name + '/log.txt','a').write("sprint statistics: %s \n\n"% (estimator.sprint_statistics()))
    
    predictions = estimator.predict(np.array(X_test,dtype=float))
    Evaluation(output_file_name,y_test,predictions,"X_test")

    pickle.dump(estimator, open(output_file_name+"/trained_automl.pkl", 'wb'))
    for i, (weight, pipeline) in enumerate(estimator.get_models_with_weights()):
        for stage_name, component in pipeline.named_steps.items():
            print(stage_name)
            if stage_name != "feature_preprocessor":
                if stage_name =='data_preprocessor':
                    # print(dir(component))
                    # print(dir(component.choice))
                    # print(component.get_preprocessor())
                    # print(component.choice.config)
                    pickle.dump(component, open(output_file_name+"/%s.pkl"%stage_name, 'wb'))
                else:
                    pickle.dump(component.choice, open(output_file_name+"/%s.pkl"%stage_name, 'wb'))
                # print(component.get_components().keys())
                
                
                open(output_file_name + '/log.txt','a').write("{0}:\n {1} : {2}\n".format(stage_name,component.choice,component.choice.get_params()))
            
    logging.info("Execution Time: %s seconds" %('{:.2f}'.format(time.time()-start_time)))    
if __name__ == '__main__':
    main()
    open(output_file_name + '/log.txt','a').write("Execution Time: %s seconds\n" %('{:.2f}'.format(time.time()-start_time)))    


import re
import json
import os
import pandas as pd

def writeData(data,writePath):
    jsondata=json.dumps(data)
    with open(writePath,'w') as f:
        json.dump(jsondata,f)

def readData(readPath):
    with open(readPath,'r') as f:
        data=json.load(f)
    data=json.loads(data)
    return data

def getVariation(data):
    
    lock=False
    index=0       
    query_name=''  
    for i in range(len(data)):
        s=data[i]
        if(s.startswith('Query=')):
            query_name=s.split(' ')[1].split('\n')[0]
        elif(s.startswith(' Score')):
            lock=True
        if(lock):
            if(s.startswith('Query')):
                index=i
                break
    lack={}
    more={}
    mut={}
    for i in range(index,len(data),4):
        if(not data[i].startswith('Query')):
            break
        q=data[i]
        ql=re.split(r' +',q)
        q_start=int(ql[1])
        q_content=ql[2].upper()

        s=data[i+2]
        sl=re.split(r' +',s)
        s_start=int(sl[1])
        s_content=sl[2].upper()
        if(len(s_content)!=len(q_content)):
            print('Error')
            break
        for j in range(len(q_content)):
            if(q_content[j]==s_content[j]):
                continue
            if(q_content[j]=='-'):
                lack[q_start+j]=s_content[j]
            elif(s_content[j]=='-'):
                more[q_start+j]=q_content[j]
            else:
                mut[q_start+j]=s_content[j]+'->'+q_content[j]
    variation={}
    variation['name']=query_name
    variation['lack']=lack
    variation['more']=more
    variation['mut']=mut
    return variation

def getMutProtein(variation,gene_data):
    
    bound=gene_data['Product_detail'].values
    gene_data['down']=0
    gene_data['up']=0  
    for i in range(len(bound)):
        down = int(re.split(r'\.+', bound[i])[0])
        up = int(re.split(r'\.+', bound[i])[-1])
        gene_data.loc[i, 'down'] = down
        gene_data.loc[i, 'up'] = up
    sequence_name=variation['name']
    subDf = gene_data[gene_data['Sequence_ID'] == sequence_name]
    if (subDf.empty):
        print('Error')
        return None
   
    lock = variation['lack']
    lock_dict = {}
    for k, v in lock.items():
        loc = int(k)
        product = subDf[subDf['down'] <= loc][subDf['up'] >= loc]['Product'].values.tolist()
        lock_dict[k] = product
    more=variation['more']
    more_dict={}
    for k,v in more.items():
        loc=int(k)
        product = subDf[subDf['down'] <= loc][subDf['up'] >= loc]['Product'].values.tolist()
        more_dict[k] = product
    mut=variation['mut']
    mut_dict={}
    for k,v in mut.items():
        loc=int(k)
        product = subDf[subDf['down'] <= loc][subDf['up'] >= loc]['Product'].values.tolist()
        mut_dict[k] = product
    mutProtein={}
    mutProtein['name']=sequence_name
    mutProtein['lock']=lock_dict
    mutProtein['more']=more_dict
    mutProtein['mut']=mut_dict
    return mutProtein

def getProteinNum(mutProtein,gene_data):
    sequence_name=mutProtein['name']
    protein_num={}
    product=gene_data[gene_data['Sequence_ID']==sequence_name]['Product'].values
    protein_num['name']=sequence_name
    for p in product:
        if(protein_num.get(p)==None):
            protein_num[p]=0
    for k,v in mutProtein['lock'].items():
        if(len(v)==0):
            print('Warning')
            continue
        for vi in v:
            if(protein_num.get(vi)==None):
                print('Waring')
                break
            protein_num[vi]+=1

    for k,v in mutProtein['more'].items():
        if(len(v)==0):
            print('Warning')
            continue
        for vi in v:
            if(protein_num.get(vi)==None):
                print('Waring')
                break
            protein_num[vi]+=1

    for k,v in mutProtein['mut'].items():
        if(len(v)==0):
            print('Warning')
            continue
        for vi in v:
            if(protein_num.get(vi)==None):
                print('Waring')
                break
            protein_num[vi]+=1
    return protein_num

def analysis(ali_path,ann_path,var_path,mut_path,pro_path):

    ali_path_list=os.listdir(ali_path)
    for file in ali_path_list:
        with open(os.path.join(ali_path,file),'r') as f:
            data=f.readlines()
        variation=getVariation(data)
        write_path=os.path.join(var_path,file.split('.')[0]+'.json')
        writeData(variation,write_path)
        ann_file_name=os.listdir(ann_path)[0]
        ann_file_path=os.path.join(ann_path,ann_file_name)
        gene_data=pd.read_csv(ann_file_path,sep='\t')
        mutProtein = getMutProtein(variation,gene_data)
        write_path = os.path.join(mut_path, file.split('.')[0] + '.json')
        writeData(mutProtein,write_path)
        proteinNum=getProteinNum(mutProtein,gene_data)
        write_path = os.path.join(pro_path, file.split('.')[0] + '.json')
        writeData(proteinNum,write_path)

def showProteinNum(pro_path):
    pro_path_list=os.listdir(pro_path)
    for file in pro_path_list:
        proteinNum=readData(os.path.join(pro_path,file))
        print(proteinNum)

def showVarNum(var_path):
    var_path_list=os.listdir(var_path)
    for file in var_path_list:
        variation=readData(os.path.join(var_path,file))
        print('lack:' + str(len(variation['lack'])), end=' ')
        print('more:' + str(len(variation['more'])), end=' ')
        print('mut:' + str(len(variation['mut'])))

def Main(types='alpha'):
    root_path=r'..\covid19_mutation'
    ali_path=os.path.join(root_path,'dataset','virus sequence alignment',types)
    ann_path=os.path.join(root_path,'dataset','virus sequence annotation',types)
    var_path = os.path.join(root_path, 'result', 'variation',types)
    mut_path=os.path.join(root_path,'result','mut2protein',types)
    pro_path=os.path.join(root_path,'result','protein2num',types)
    if(not os.path.exists(var_path)):
        os.mkdir(var_path)
    if(not os.path.exists(mut_path)):
        os.mkdir(mut_path)
    if(not os.path.exists(pro_path)):
        os.mkdir(pro_path)

    analysis(ali_path,ann_path,var_path,mut_path,pro_path)
    showVarNum(var_path)
    showProteinNum(pro_path)

covids=['alpha','beta','gamma','delta','omicron']
for covid in covids:
    Main(covid)

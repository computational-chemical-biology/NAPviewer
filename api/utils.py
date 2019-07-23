import pandas as pd 
import json
import numpy
import requests
import colorsys
import re

# https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
def get_N_HexCol(N=5):
    HSV_tuples = [(x * 1.0 / N, 0.5, 0.5) for x in range(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x * 255), colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#%02x%02x%02x' % tuple(rgb))
    return hex_out


def downloadRes(jobid):
    fname = jobid+'.zip'
    url = "https://proteomics2.ucsd.edu/ProteoSAFe/DownloadResult?task=" + jobid + "&view=summary_report"

    response = requests.post(url)

    with open(fname, "wb") as output_file:
        for block in response.iter_content(1024):
            output_file.write(block)

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

# add connected component id and name
def formattable(jobid):
    tab = pd.read_table('api/static/downloads/'+jobid+'/node_attributes_table.tsv') 
    #tab = tab.tail()
    tab2 = tab[['cluster.index', 'parent.mass', 'RTMean', 'LibraryID', 'FusionID', 'ConsensusID', 'MetFragID', 'MetFragScore' ]] 
    tab2 = tab2.rename(columns = {'cluster.index': 'clusterindex', 'parent.mass': 'parentmass', 'MetFragScore':'CompleteRes'}) 
    tab2.fillna({'FusionID': '', 'ConsensusID': '', 'MetFragID': '', 'CompleteRes': ''}, inplace=True)
    tab2.reset_index(inplace=True, drop=True)  

    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_annotations_DB#%7B%22main.Compound_Name_input%22%3A%22'
    eurl = '%22%7D'
    gid = tab.loc[0,'ProteoSAFeClusterLink'].split('&')[0].split('=')[1]

    rjobid = jobid[:5] 

    for i in range(tab2.shape[0]):
        if type(tab2['LibraryID'][i]) == float:
            tab2.loc[i,'LibraryID'] = '<a href=\"'+ tab.loc[i,'ProteoSAFeClusterLink'] + '\" target=\"_blank\">'+'N/A'+'</a>'
        else:	    
            tab2.loc[i,'LibraryID'] = '<a href=\"'+ burl+gid+purl+str(tab2.loc[i,'LibraryID'])+eurl+ '\" target=\"_blank\">'+str(tab2.loc[i,'LibraryID'])+'</a>'

        if tab2['MetFragID'][i] != '':
            tab2.loc[i,'MetFragID'] = tab2['MetFragID'][i].split(',')[0]   
            tab2.loc[i,'MetFragID'] = ('<a href=\"javascript:showColumn2(\'G\', \'M\', \'' +  
                                        tab2['MetFragID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+
					');\">G</a>'+'\n'+
                                        '<a href=\"javascript:showColumn2(\'S\', \'M\', \'' +  
                                        tab2['MetFragID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+
					');\">S</a>'+'\n'+
                                        '<a href=\"javascript:showColumn2(\'F\', \'M\', \'' +  
                                        tab2['MetFragID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+

					');\">F</a>')
        if tab2['FusionID'][i] != '':
            tab2.loc[i,'FusionID'] = tab2['FusionID'][i].split(',')[0]   
            tab2.loc[i,'FusionID'] = ('<a href=\"javascript:showColumn1(\'G\', \'F\', \'' +  
                                        tab2['FusionID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+
					');\">G</a>'+'\n'+
                                        '<a href=\"javascript:showColumn1(\'S\', \'F\', \'' +  
                                        tab2['FusionID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+
					');\">S</a>'+'\n'+
                                        '<a href=\"javascript:showColumn1(\'F\', \'F\', \'' +  
                                        tab2['FusionID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+

					');\">F</a>')

        if tab2['ConsensusID'][i] != '':
            tab2.loc[i,'ConsensusID'] = tab2['ConsensusID'][i].split(',')[0]   
            tab2.loc[i,'ConsensusID'] = ('<a href=\"javascript:showColumn1(\'G\', \'C\', \'' +  
                                        tab2['ConsensusID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+
					');\">G</a>'+'\n'+
                                        '<a href=\"javascript:showColumn1(\'S\', \'C\', \'' +  
                                        tab2['ConsensusID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+
					');\">S</a>'+'\n'+
                                        '<a href=\"javascript:showColumn1(\'F\', \'C\', \'' +  
                                        tab2['ConsensusID'][i] + '\', \''+ str(tab2['clusterindex'][i]) + '\','+
                                        '\''+ rjobid + '\''+

					');\">F</a>')

   
        if tab2['CompleteRes'][i] != '':

            url =  'javascript:showColumn1(\'F\', \'L\',' + '\'all\'' + ', \''+ str(tab2['clusterindex'][i]) + '\', \''+ rjobid + '\');'
            tab2.loc[i,'CompleteRes'] = '<a href=\"'+ url + '\">Link</a>'

    meas = [tab2.iloc[i].to_dict() for i in range(tab2.shape[0])] 
    meas = {'data': meas}
    return meas 


def linkpattern(id):

    if bool(re.match('^CCMSLIB', id)): 
        url = "http://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^HMDB', id)):
        url = "http://www.hmdb.ca/metabolites/"+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^SN', id)):
        url = "http://bioinf-applied.charite.de/supernatural_new/index.php?site=search5&np_id="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^CHEBI', id)):
        url = "https://www.ebi.ac.uk/chebi/searchId.do?chebiId="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^CID', id)):
        url = "https://pubchem.ncbi.nlm.nih.gov/compound/"+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^Chemical-Structure\\.', id)):
        url = "http://www.chemspider.com/"+id+".html"
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^cs', id)):
    #elif bool(re.match('^\d{3,}', id)):
        # provisory for marinlit
        #id2 = id.zfill(12)
        #url = "http://pubs.rsc.org/marinlit/compound/cs"+id2
        url = "http://pubs.rsc.org/marinlit/compound/"+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^FDB', id)):
        url = 'http://foodb.ca/compounds/'+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^DB', id)):
        url = 'https://www.drugbank.ca/drugs/'+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^BGC\d{7,}', id)):
        url = 'https://mibig.secondarymetabolites.org/repository/'+id+'/index.html#cluster-1'
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^NOR\d{5,}', id)):
        url = 'http://bioinfo.lifl.fr/norine/result.jsp?ID='+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('[A-Z]{3,}', id)):
        # provisory for marinlit
        url = 'http://dnp.chemnetbase.com/faces/chemical/FullScreenEntry.xhtml?product=DNP&fromExport=falsePOSI'+id+'&id='+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    else:
        #return "http://dnp.chemnetbase.com/AAA00.entry?parentCHNumber="+id+"&exno="+id
        # http://dnp.chemnetbase.com/faces/chemical/FullScreenEntry.xhtml?product=DNP&fromExport=falsePOSIBBB04&id=BBB04
        url = "https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds&query_type=text&query="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'

#    '<a href=\"javascript:showColumn1(\' world\'); 
#     javascript:showColumn2(\' world\');\" >'+tab2['MetFragID'][0]+ '</a>'
def formatsubtable(jobid, clustindex, mlist, tabgnps, lid):
    fname = 'api/static/downloads/'+jobid
    #with open(fname+'/mlist.json', 'r') as f:
    #    mlist = json.load(f)

    #tabgnps = pd.read_table(fname+'/tabgnps.tsv', sep='\t')
    pos = numpy.where(tabgnps['cluster.index']==int(clustindex))[0][0]
    mtab = pd.DataFrame(mlist[pos])
    # recover instances where we 
    # have fusion/consensus
    if mtab.shape[0] < 4:
        mtab = pd.DataFrame(lid[pos])[['Identifier', 'MonoisotopicMass', 'NoExplPeaks', 'Score']]
	
    mtab = mtab.rename(columns = {'fusion2': 'Consensus', 'fusion': 'Fusion'}) 

    if not 'Consensus' in mtab.columns:
        mtab['Consensus'] = ''
    if not 'Fusion' in mtab.columns:
        mtab['Fusion'] = ''
    if not 'MCSS' in mtab.columns:
        mtab['MCSS'] = ''
    if not 'superclass_name' in mtab.columns:
        mtab['superclass_name'] = ''
    if not 'class_name' in mtab.columns:
        mtab['class_name'] = ''

    lcol = list(set(mtab['MCSS'])) 
    if len(lcol)>1:
        col = get_N_HexCol(len(lcol))
        dcol = dict(zip(lcol, col))	    

    rjobid = jobid[:5] 
    for i in range(mtab.shape[0]):
        if mtab['NoExplPeaks'][i] != 0:
            mtab.loc[i,'NoExplPeaks'] = ('<a href=\"javascript:showColumn2(\'F\', \'M\',\'' +  
                                        re.sub('\s+', '', mtab['Identifier'][i]) + '\', \''+ str(clustindex) + '\','+

                                        '\''+ rjobid + '\''+
					');\">'+ mtab.loc[i,'NoExplPeaks'] +'</a>'
                                      )
    for j in range(mtab.shape[0]):
        mtab.loc[j,'Identifier'] = linkpattern(str(mtab.loc[j,'Identifier']))
        if str(mtab['MCSS'][j])!='':
            mtab.loc[j, 'MCSS'] = ('<span style=\"background-color:'+dcol[mtab.loc[j, 'MCSS']]+
	                           '\"><font color="white">G'+str(mtab.loc[j, 'MCSS'])+'</font></span>'
				   )
    
    meas = [mtab.iloc[i].to_dict() for i in range(mtab.shape[0])] 
    meas = {'data': meas}
    return meas 





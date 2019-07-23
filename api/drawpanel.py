
import re, os, json
from rdkit import Chem
from rdkit.Chem import Draw
from bokeh.plotting import figure, output_file, save, ColumnDataSource
from bokeh.models import HoverTool
from pyteomics import mgf
import pandas as pd
import numpy as np
from api import visJS_module
#from visJS2jupyter import visJS_module
import networkx as nx
import unicodedata

# https://stackoverflow.com/questions/517923/what-is-the-best-way-to-remove-accents-in-a-python-unicode-string
def remove_accents(input_str):
    nfkd_form = unicodedata.normalize('NFKD', input_str)
    return u"".join([c for c in nfkd_form if not unicodedata.combining(c)])

def plotPannel(tabgnps, idtab, clusterid, score, dr, nstruct=10):
#tabgnps = pd.read_table('tabgnps.tsv', sep='\t')
#tabgnps.head()
#with open('lid.json', 'r') as f:
#idtab = json.load(f)
    """Plots a pannel with candidate structures.

    Retrieves a list of candidate structures
    from each method and order by its score.

    Args:
        tabgnps: Node attributes table, with node order (do we need that?).
        idtab: a dict list with all node candidades (should we import the whole list?).
            to fetch.
        clusterid: cluster id (row index?) to select the dict.
        score: Scoring method to order the table.
        nstruct: Number of structures to display.

    Returns:
	A png pannel with the structures ordered by score, with structure id and score in
	the plot.
        example:

    Raises:
        IOError: ...
    """
#clusterid = 5559 
    filename = dr + str(clusterid)+score+'.png'

    if(os.path.exists(filename)):
        pass

    clsid = np.where(tabgnps['cluster.index']==int(clusterid))[0][0]
    attr = pd.DataFrame(idtab[clsid])
    attr[score] = attr[score].apply(pd.to_numeric)
    attr.sort_values(by=score, ascending=False, inplace=True)
    attr.reset_index(inplace=True)

    if attr.shape[0] > nstruct:
        attr = attr.head(nstruct)

    mols = [Chem.MolFromSmiles(x) for x in list(attr["SMILES"])]
    leg = []
#score = 'Score'
    for i in range(attr.shape[0]):
        cscore = str(round(float(attr[score][i]), 3))
        leg.append(remove_accents(attr["Identifier"][i])+'-'+cscore)

    img=Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(200,200),legends=leg)
    img.save(filename)
    return '<html><img src=\"'+filename+'\" style="width:100%;height:100%;"></html>'

#scanindex = 5805 
#mgffile = 'allspectra.mgf' 
#create a dictionary with only spectra 
# in the tabgnps
def readSpectrum(mgffile, scanindex):
    msms1 = []
    with mgf.read(mgffile) as allspectra:
        for spectrum in allspectra:
            n = int(re.sub("\D+", "", spectrum['params']['title'] ))
            if n == scanindex:
                mz = spectrum['m/z array']
                inty = spectrum['intensity array']
                for i in range(len(mz)):
                    msms1.append((mz[i], inty[i]))
    return(pd.DataFrame(msms1, columns=['Mass', 'Intensity']))

def getSpectra(mgffile, tabgnps):
    spectra = {}
    #tabgnps['cluster.index'] = tabgnps['cluster.index'].astype(int)
    with mgf.read(mgffile) as allspectra:
        for spectrum in allspectra:
            msms1 = []
            n = int(re.sub("\D+", "", spectrum['params']['title'] ))
            idx = np.where(tabgnps['cluster.index']==n)[0]
            if len(idx):
                mz = spectrum['m/z array']
                inty = spectrum['intensity array']
                for i in range(len(mz)):
                    msms1.append((mz[i], inty[i]))
                spectra[tabgnps['cluster.index'][idx[0]]] = pd.DataFrame(msms1, columns=['Mass', 'Intensity'])
    return(spectra)

def getSingleSpectrum(mgffile, scanindex):
    spectrum = {}
    f = open(mgffile)
    lines = np.array(f.readlines())
    f.close()
    bg = np.where(lines=='BEGIN IONS\n')
    ed = np.where(lines=='END IONS\n')
    p = np.where(lines=='SCANS='+str(scanindex)+'\n')
    bgp = bg[0][np.where((bg[0]-p[0])>0)[0][0]-1]
    edp = ed[0][np.where((ed[0]-p[0])>0)[0][0]]+1
    msms1 = []
    for line in lines[bgp:edp]:
        try:
            mz, rt = line.split(' ')
            msms1.append((float(mz), float(rt)))
        except:
            pass
    spectrum[scanindex] = pd.DataFrame(msms1, columns=['Mass', 'Intensity'])
    return(spectrum)

def drawSingleMol(smi, fout):
    mol = Chem.MolFromSmiles(str(smi))
    Draw.MolToFile(mol, fout)

#cpdid = 26624 
# Fix no ExplPeaks exception
def plotFragSpectrum(mgffile, tabgnps, insilicopred, scanindex, cpdid, fragdir):
    pos = np.where(tabgnps['cluster.index']==scanindex)[0][0]
    insilicopred = insilicopred[pos]
    insilicopred = pd.DataFrame(insilicopred)
    # Do we have identifiers with space in the names?
    insilicopred['Identifier'] = [re.sub('\s+', '', x) for x in insilicopred['Identifier']]
    formulasexplpeaks = insilicopred[insilicopred['Identifier'].astype(str)==str(cpdid)]['FormulasOfExplPeaks'].iloc[0].split(';')
    formulasexplpeaks = pd.DataFrame([x.split(':') for x in formulasexplpeaks])
    smilesofexplpeaks = insilicopred[insilicopred['Identifier'].astype(str)==str(cpdid)]['SmilesOfExplPeaks'].iloc[0].split(';')
    smilesofexplpeaks = pd.DataFrame([x.split(':')[1] for x in smilesofexplpeaks])
    #spectrum = readSpectrum(mgffile, scanindex) 
    spectrum = mgffile[scanindex]
    specinfo = pd.concat([formulasexplpeaks, smilesofexplpeaks], axis=1, ignore_index=True)
    specinfo.columns = ['Massref', 'Formula', 'Frag']
    specinfo['Massref'] = specinfo['Massref'].astype(float)
    specinfo = pd.merge(spectrum, specinfo, left_on='Mass', right_on='Massref', how='left')
    specinfo.fillna('', inplace=True)

    # folder for all fragments, and png with scanindex
    #fragdir = str(scanindex) + '_' + str(cpdid) + '_frags' 
    specinfo['Link'] = ''
    specinfo['Lab'] = ''
    for x in range(specinfo.shape[0]):
        if specinfo['Frag'][x]!='':
            try:
                fign = str(scanindex) + '_' + str(cpdid) + '_' + str(x)+'.png'
                drawSingleMol(specinfo['Frag'][x], fragdir+'/'+fign)
                specinfo.loc[x, 'Link'] = fign
#specinfo.loc[x, 'Lab'] = '\n'+str(round(specinfo['Mass'][x], 3))+'\n'+specinfo['Formula'][x]
                specinfo.loc[x, 'Lab'] = specinfo['Formula'][x]
            except:
                pass

#specinfo['Link'] = specinfo['Link'].astype(str) 
    idx = specinfo["Link"]!=''
    specinfo2 = specinfo[idx]

    source = ColumnDataSource(
        data=dict(
		x=list(specinfo2["Mass"]),
		y=list(specinfo2["Intensity"]),
		desc=list(specinfo2["Lab"]),
		imgs = list(specinfo2["Link"]),
        )
	)

    p = figure(sizing_mode='stretch_both', tools="pan,wheel_zoom,box_zoom,reset,save",
              title="Mouse over the dots")

    x = list(specinfo["Mass"])
    y = list(specinfo["Intensity"])

    line = p.multi_line([[i,i] for i in x], [[0,i] for i in y] )
    circle = p.circle('x', 'y', size=5, source=source)

    hover = HoverTool(
	    tooltips="""
	    <div>
		<div>
		    <img
			src="@imgs" height="150" alt="@imgs" width="150"
			style="float: left; margin: 0px 15px 15px 0px;"
			border="2"
		    ></img>
		</div>
		<div>
		    <span style="font-size: 17px; font-weight: bold;">@desc</span>
		    <span style="font-size: 15px; color: #966;">[$index]</span>
		</div>
		</div>
	    </div>
	    """,
	    renderers=[circle]
	)

    p.add_tools(hover)

    idx2 = [x for x in range(specinfo.shape[0]) if specinfo["Lab"][x]!='' and  specinfo["Link"][x]=='']

    if len(idx2):
        specinfo2 = specinfo.iloc[idx2]
        source2 = ColumnDataSource(
	    data=dict(
		x=list(specinfo2["Mass"]),
		y=list(specinfo2["RT"]),
		desc=list(specinfo2["Lab"]),
		imgs = list(specinfo2["Link"]),
	    )
	)
        circle2 = p.circle('x', 'y', size=5, source=source2)
        hover2 = HoverTool(
	    tooltips="""
	    <div>
		 <div>
		    <span style="font-size: 17px; font-weight: bold;">@desc</span>
		    <span style="font-size: 15px; color: #966;">[$index]</span>
		</div>
		</div>
	    </div>
	    """,
	    renderers=[circle2]
	)
        p.add_tools(hover2)

    save(p, fragdir+str(scanindex) + '_' + str(cpdid)+'.html')

#net = pd.read_table('net.tsv', sep='\t')
#scanindex = 3448 
#method = 'MF'

def plotGraph(tabgnps, idtab, net, scanindex, method, dr, option=2):
    # option 1: connected component
    # option 2: direct neighbors
    if option==1:
        c1 = net['V1']==scanindex
        c2 = net['V2']==scanindex
        id = np.where([a or b for a, b in zip(c1, c2)])[0][0]
        subnet = net[net['V7']==net.iloc[id]['V7']]
    elif option==2:
        c1 = net['V1']==scanindex
        c2 = net['V2']==scanindex
        id = list(np.where([a or b for a, b in zip(c1, c2)])[0])
        subnet = net.iloc[id]

    G = nx.from_pandas_dataframe(subnet,source='V1',target='V2',edge_attr = 'V5')

    nodes = list(G.nodes())
    edges = list(G.edges())
    clpos = [np.where(tabgnps['cluster.index']==x)[0][0] for x in nodes]

    image = []
    for c in clpos:
        if type(idtab[c]) is not list:
            image.append('')
            continue

        if type(tabgnps['Smiles'][c]) is str:
            smi = tabgnps['Smiles'][c]
            img = 'LM'+str(tabgnps['cluster.index'][c])+'.png'
            fout = dr+img
        else:
            tmpidtab = pd.DataFrame(idtab[c])
            if method=='MF' or (('fusion' not in tmpidtab.columns ) and ('fusion2' not in tmpidtab.columns )):
                smi = tmpidtab.iloc[0]['SMILES']
                img = 'MF'+str(tabgnps['cluster.index'][c])+'.png'
            elif (method=='FS' and ('fusion' in tmpidtab.columns)) or (method=='CS' and ('fusion2' not in tmpidtab.columns)):
                tmpidtab.sort_values(by='fusion', ascending=False, inplace=True)
                smi = tmpidtab.iloc[0]['SMILES']
                img = 'FS'+str(tabgnps['cluster.index'][c])+'.png'
            elif method=='CS':
                tmpidtab.sort_values(by='fusion2', ascending=False, inplace=True)
                smi = tmpidtab.iloc[0]['SMILES']
                img = 'CS'+str(tabgnps['cluster.index'][c])+'.png'

            fout = dr+img
	# check the figures in the folder
        try:
            drawSingleMol(smi=smi, fout=fout)
            image.append(img)
        except:
            image.append('')

    pos = nx.spring_layout(G)

    nodes_dict = []

    for j in range(len(nodes)):
        n = nodes[j]
        label = ('id:'+str(tabgnps['cluster.index'][clpos[j]])+' '+
            'm/z:'+str(round(tabgnps['parent.mass'][clpos[j]], 3)))
        if image[j]!='':
            shape = 'circularImage'
            address = image[j]
            ibool = True
            if 'MF' in address:
                nbcol = '#00FFFF'
            elif 'FS' in address:
                nbcol = '#0000FF'
            elif 'CS' in address:
                nbcol = '#FF0000'
            else:
                nbcol = '#00FF00'
            nbwidth = 20
        else:
            shape = 'triangle'
            address = 'undefined'
            ibool = False
            nbcol = '#FF0000'
            nbwidth = 1
        nodes_dict.append(
          {"id":n,
           #"degree":nx.degree(G,n),
           "degree": 1,
           "border_width": nbwidth,
           "label":label,
           "title":label,
           "node_size": 100,
           "node_shape_use_border_with_image":ibool,
           "border_color": nbcol,
           "node_shape": shape,
            "node_image":address,
            "x":pos[n][0]*200,
            "y":pos[n][1]*200
           }
        )

    node_map = dict(zip(nodes,range(len(nodes))))  # map to indices for source/target in edges

    # still can't change edge label
    elab = str(0.75)
    edges_dict = [{"source":node_map[edges[i][0]], "target":node_map[edges[i][1]],
                 "color":"gray","edge_label":0.75, "edge_width": 5} for i in range(len(edges))]

    f = open(dr+str(scanindex)+method+'.html', 'w+')
    f.write(
    # set some network-wide styles
    visJS_module.visjs_network(nodes_dict,edges_dict,
                   node_size_multiplier=50,
                   #node_color_border= '#0000FF',
                   #graph_width = '100%',
                   #graph_height = '100%',
                   node_size_transform = 'Math.abs',
                   node_font_size = 50,
                   edge_width=5,
                   navigation_buttons = True,
                   showButton = True,
                   output='html')
    )
    f.close()




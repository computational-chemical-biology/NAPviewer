#!/home/flaskappuser/miniconda2/envs/napviewer/bin/python
from flask import Flask, render_template, request, redirect, url_for, make_response
#from flask_session import Session

import pandas as pd
import json
import numpy as np
from api.utils import formattable, MyEncoder, downloadRes, formatsubtable
import os, re

from api.drawpanel import *

app = Flask(__name__)

app.config['SECRET_KEY'] = 'super secret key'
#app.config['DEBUG'] = True

# https://stackoverflow.com/questions/18967441/add-a-prefix-to-all-flask-routes/36033627#36033627
class PrefixMiddleware(object):

    def __init__(self, app, prefix=''):
        self.app = app
        self.prefix = prefix

    def __call__(self, environ, start_response):

        if environ['PATH_INFO'].startswith(self.prefix):
            environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
            environ['SCRIPT_NAME'] = self.prefix
            return self.app(environ, start_response)
        else:
            start_response('404', [('Content-Type', 'text/plain')])
            return ["This url does not belong to the app.".encode()]




#SESSION_TYPE = 'filesystem'
#app.config.from_object(__name__)
#Session(app)

app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix='/NAPviewer')

@app.route('/')
def page_view():
    jobid = str(request.args.get('task'))
    fname = 'api/static/downloads/'+jobid
    if not os.path.exists(fname):
        return redirect(url_for('download', jobid=jobid))
    else:
        return render_template('view.html', jobid=jobid)

@app.route('/download')
def download():
    jobid = str(request.args.get('jobid'))
    fname = 'api/static/downloads/'+jobid

    if not os.path.exists(fname):
        #downloadRes(jobid)
        #os.system('Rscript rda2json.R '+jobid+'.zip '+fname)
        os.system('Rscript api/rda2json2.R '+jobid+' '+fname)
        fname = 'api/static/downloads/'+jobid+'/panel/'
        os.mkdir(fname)
        fname = 'api/static/downloads/'+jobid+'/nodeimages/'
        os.mkdir(fname)
        fname = 'api/static/downloads/'+jobid+'/frags/'
        os.mkdir(fname)
        return redirect(url_for('page_view', task=jobid))
    else:
        return redirect(url_for('page_view', task=jobid))


@app.route('/loadtable')
def load():
    # display name of first candidate

    jobid = str(request.args.get('jobid'))

    meas = formattable(jobid)
    #return jsonify({'data': meas}) 
    return json.dumps(meas, cls=MyEncoder)

@app.route('/loadsubtable')
def subload():
    clsid = str(request.args.get('clsid'))
    jobid = str(request.args.get('jobid'))

    fname = 'api/static/downloads/'+jobid
    tabgnps = pd.read_table(fname+'/tabgnps.tsv', sep='\t')
    with open(fname+'/lid.json', 'r') as f:
        lid = json.load(f)
    with open(fname+'/mlist.json', 'r') as f4:
        mlist = json.load(f4)

    meas = formatsubtable(jobid, clsid, mlist, tabgnps, lid)
    #return jsonify({'data': meas}) 
    return json.dumps(meas, cls=MyEncoder)


@app.route("/plot2")
def plot2():
    tpid = str(request.args.get('type'))
    score = str(request.args.get('score'))
    strid = str(request.args.get('strid'))
    clsid = str(request.args.get('clsid'))
    rjobid = str(request.args.get('rjobid'))

    fls = os.listdir('api/static/downloads')
    jobid = [x for x in fls if bool(re.search(rjobid, x))]

    if len(jobid)==1:
        jobid = jobid[0]
    else:
        hd = '<h3>Duplicated taskid, please contact the maintainer</h3>'
        return render_template('plot.html', image='<img src=\"\">',
                               header=hd)

    if tpid == 'S':
        if score == 'M':
            fname = 'static/downloads/'+jobid+'/panel/'
            filename = fname + str(clsid)+'Score'+'.png'
            if os.path.exists(filename):
                hd = '<h3>Node '+ str(clsid) + ' MetFrag candidates</h3>'
                return render_template('plot.html', image='<img src=\"'+filename+'\" >',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                with open(fn+'/lid.json', 'r') as f:
                    lid = json.load(f)

                fname = 'api/static/downloads/'+jobid+'/panel/'
                plotPannel(tabgnps, lid, clsid, score='Score', dr=fname, nstruct=10)
                fname = 'static/downloads/'+jobid+'/panel/'
                filename = fname + str(clsid)+'Score'+'.png'
                hd = '<h3>Node '+ str(clsid) + ' MetFrag candidates</h3>'
                return render_template('plot.html', image='<img src=\"'+filename+'\">',
                       header=hd)
        elif score == 'F':
            fname = 'static/downloads/'+jobid+'/panel/'
            filename = fname + str(clsid)+'fusion'+'.png'
            if os.path.exists(filename):
                hd = '<h3>Node '+ str(clsid) + ' Fusion candidates</h3>'
                return render_template('plot.html', image='<img src=\"'+filename+'\" >',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                with open(fn+'/fusion.json', 'r') as f2:
                    fus = json.load(f2)
                fname = 'api/static/downloads/'+jobid+'/panel/'
                plotPannel(tabgnps, fus, clsid, score='fusion', dr=fname, nstruct=10)
                fname = 'static/downloads/'+jobid+'/panel/'
                filename = fname + str(clsid)+'fusion'+'.png'
                hd = '<h3>Node '+ str(clsid) + ' Fusion candidates</h3>'
                return render_template('plot.html', image='<img src=\"'+filename+'\">',
                       header=hd)

        elif score == 'C':
            fname = 'static/downloads/'+jobid+'/panel/'
            filename = fname + str(clsid)+'fusion2'+'.png'
            if os.path.exists(filename):
                hd = '<h3>Node '+ str(clsid) + ' Consensus candidates</h3>'
                return render_template('plot.html', image='<img src=\"'+filename+'\" >',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                with open(fn+'/consensus.json', 'r') as f3:
                    cons = json.load(f3)
                fname = 'api/static/downloads/'+jobid+'/panel/'
                plotPannel(tabgnps, cons, clsid, score='fusion2', dr=fname, nstruct=10)
                fname = 'static/downloads/'+jobid+'/panel/'
                filename = fname + str(clsid)+'fusion2'+'.png'
                hd = '<h3>Node '+ str(clsid) + ' Consensus candidates</h3>'
                return render_template('plot.html', image='<img src=\"'+filename+'\">',
                       header=hd)

    elif tpid == 'G':
        if score == 'M':
            fname = 'api/static/downloads/'+jobid+'/nodeimages/'
            fle = fname+str(clsid)+'MF'+'.html'
            if os.path.exists(fle):
                fname = 'static/downloads/'+jobid+'/nodeimages/'
                fle = fname+str(clsid)+'MF'+'.html'
                hd = '<h3>Node '+ str(clsid) + ' MetFrag direct neighbors</h3>'
                return render_template('plot2.html', image='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                net = pd.read_table(fn+'/net.tsv', sep='\t')
                with open(fn+'/lid.json', 'r') as f:
                    lid = json.load(f)
                plotGraph(tabgnps, lid, net, int(clsid), method='MF', dr=fname, option=2)
                fname = 'static/downloads/'+jobid+'/nodeimages/'
                fle = fname+str(clsid)+'MF'+'.html'
                hd = '<h3>Node '+ str(clsid) + ' MetFrag direct neighbors</h3>'
                return render_template('plot2.html',
                       image='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)
        elif score == 'F':
            fname = 'api/static/downloads/'+jobid+'/nodeimages/'
            fle = fname+str(clsid)+'FS'+'.html'
            if os.path.exists(fle):
                fname = 'static/downloads/'+jobid+'/nodeimages/'
                fle = fname+str(clsid)+'FS'+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Fusion direct neighbors</h3>'
                return render_template('plot2.html', image='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                net = pd.read_table(fn+'/net.tsv', sep='\t')
                with open(fn+'/fusion.json', 'r') as f2:
                    fus = json.load(f2)
                plotGraph(tabgnps, fus, net, int(clsid), method='FS', dr=fname, option=2)
                fname = 'static/downloads/'+jobid+'/nodeimages/'
                fle = fname+str(clsid)+'FS'+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Fusion direct neighbors</h3>'
                return render_template('plot2.html',
                       image='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                      header=hd)
        elif score == 'C':
            fname = 'api/static/downloads/'+jobid+'/nodeimages/'
            fle = fname+str(clsid)+'CS'+'.html'
            if os.path.exists(fle):
                fname = 'static/downloads/'+jobid+'/nodeimages/'
                fle = fname+str(clsid)+'CS'+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Consensus direct neighbors</h3>'
                return render_template('plot2.html', image='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                net = pd.read_table(fn+'/net.tsv', sep='\t')
                with open(fn+'/consensus.json', 'r') as f3:
                    cons = json.load(f3)
                plotGraph(tabgnps, cons, net, int(clsid), method='CS', dr=fname, option=2)
                fname = 'static/downloads/'+jobid+'/nodeimages/'
                fle = fname+str(clsid)+'CS'+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Consensus direct neighbors</h3>'
                return render_template('plot2.html',
                       image='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)

    elif tpid == 'F':
        if score == 'M':
            fname = 'api/static/downloads/'+jobid+'/frags/'
            fle = fname+str(clsid)+'_'+str(strid)+'.html'
            if os.path.exists(fle):
                fname = 'static/downloads/'+jobid+'/frags/'
                fle = fname+str(clsid)+'_'+str(strid)+'.html'
                hd = '<h3>Node '+ str(clsid) + ' MetFrag ' + str(strid) + ' pred. fragments</h3>'
                return render_template('plot3.html', frag='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                #allspec = getSpectra(fn+'/allspectra.mgf', tabgnps)
                allspec = getSingleSpectrum(fn+'/allspectra.mgf', int(clsid))
                with open(fn+'/lid.json', 'r') as f:
                    lid = json.load(f)
                fname = 'api/static/downloads/'+jobid+'/frags/'
                plotFragSpectrum(allspec, tabgnps, lid, int(clsid), strid, fname)
                fname = 'static/downloads/'+jobid+'/frags/'
                fle = fname+str(clsid)+'_'+str(strid)+'.html'
                hd = '<h3>Node '+ str(clsid) + ' MetFrag ' + str(strid) + ' pred. fragments</h3>'
                return render_template('plot3.html',
                        frag='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                        header=hd)
        # Only the header is changing here
        # if the call comes from cand. table, add a parameter to change
        # the header
        elif score == 'F':
            fname = 'api/static/downloads/'+jobid+'/frags/'
            fle = fname+str(clsid)+'_'+str(strid)+'.html'
            if os.path.exists(fle):
                fname = 'static/downloads/'+jobid+'/frags/'
                fle = fname+str(clsid)+'_'+str(strid)+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Fusion ' + str(strid) + ' pred. fragments</h3>'
                return render_template('plot3.html', frag='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                #allspec = getSpectra(fn+'/allspectra.mgf', tabgnps)
                allspec = getSingleSpectrum(fn+'/allspectra.mgf', int(clsid))
                with open(fn+'/lid.json', 'r') as f:
                    lid = json.load(f)
                fname = 'api/static/downloads/'+jobid+'/frags/'
                plotFragSpectrum(allspec, tabgnps, lid, int(clsid), strid, fname)
                fname = 'static/downloads/'+jobid+'/frags/'
                fle = fname+str(clsid)+'_'+str(strid)+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Fusion ' + str(strid) + ' pred. fragments</h3>'
                return render_template('plot3.html',
                       frag='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)
        elif score == 'C':
            fname = 'api/static/downloads/'+jobid+'/frags/'
            fle = fname+str(clsid)+'_'+str(strid)+'.html'
            if os.path.exists(fle):
                fname = 'static/downloads/'+jobid+'/frags/'
                fle = fname+str(clsid)+'_'+str(strid)+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Consensus ' + str(strid) + ' pred. fragments</h3>'
                return render_template('plot3.html', frag='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)
            else:
                fn = 'api/static/downloads/'+jobid
                tabgnps = pd.read_table(fn+'/tabgnps.tsv', sep='\t')
                #allspec = getSpectra(fn+'/allspectra.mgf', tabgnps)
                allspec = getSingleSpectrum(fn+'/allspectra.mgf', int(clsid))
                with open(fn+'/lid.json', 'r') as f:
                    lid = json.load(f)
                fname = 'api/static/downloads/'+jobid+'/frags/'
                plotFragSpectrum(allspec, tabgnps, lid, int(clsid), strid, fname)
                fname = 'static/downloads/'+jobid+'/frags/'
                fle = fname+str(clsid)+'_'+str(strid)+'.html'
                hd = '<h3>Node '+ str(clsid) + ' Consensus ' + str(strid) + ' pred. fragments</h3>'
                return render_template('plot3.html',
                       frag='<object width="100%" height="100%" data=\"'+fle+'\"></object>',
                       header=hd)

    if strid=='all':
        return render_template('plot4.html', clustid=clsid, jobid=jobid)


@app.route("/getimage")
def getimage():
    path = str(request.args.get('path'))
    #resp = make_response(open(path.encode('utf-8').strip()).read())
    #resp.content_type = "image/jpeg"
    fullpath = 'api/static/'+jobid+'/nodeimages/'+path
    return app.send_static_file(path)

#sess = Session()
#sess.init_app(app)

if __name__ == '__main__':
    #app.run(host='137.110.133.15', port=5001)
    #app.run(host='127.0.0.1')
    app.run(port=5001)


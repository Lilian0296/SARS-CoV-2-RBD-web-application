from flask import Flask, render_template, request, session, redirect, url_for, flash
import pandas as pd
import string as s
import random
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import plotly.express as px
from kaleido.scopes.plotly import PlotlyScope
import plotly
import plotly.express as px
import plotly.graph_objects as go
import chart_studio
import plotly_utils
import re
import pygal
import os
import json
from werkzeug.utils import secure_filename

pd.set_option('display.max_colwidth', 40)
# index for matrix
c_index=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
c_number=list(range(0,20))
r_positions=list(range(331,532))
r_number=list(range(0,201))

col=dict(zip(c_index,c_number))
row=dict(zip(r_number,r_positions))

    
# matrix by wuhan fill na with median
data=pd.read_csv("data_fillnamedian.csv")
data=data[data["target"]=="Wuhan-Hu-1"]
wuhan_m=[]
for i in range(len(r_positions)):
    item=np.array(data[data["position"]==r_positions[i]]["bind"])
    wuhan_m.append(item)
wuhan_m=np.array(wuhan_m)    



# Define folder to save uploaded files to process further
UPLOAD_FOLDER = os.path.join('static', 'uploads')

# Define allowed files (for this example I want only csv file)
ALLOWED_EXTENSIONS = {'fasta'}
 
 
app = Flask(__name__)

# Configure upload file path flask
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

#It will allow below 16MB contents only, you can change it
app.config['MAX_CONTENT_LENGTH'] = 2048 * 1024 * 1024
 
# Define secret key to enable session
app.secret_key = 'This is your secret key to utilize session in Flask'

@app.route('/')
def home():
    return render_template('RBD_home.html')



@app.route('/',  methods=("POST", "GET"))
def uploadFile():
    if request.method == 'POST':
        if 'uploaded-file[]' not in request.files:
            flash('No file part')
            return redirect(request.url)
        # upload file flask
        uploaded_df = request.files.getlist('uploaded-file[]')
        
        for file in uploaded_df:
             # Extracting uploaded data file name
            data_filename = secure_filename(file.filename)
            # flask upload file to database (defined uploaded folder in static path)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], data_filename))
            # Storing uploaded file path in flask session
        flash('File(s) successfully uploaded')
        return render_template('RBD_home.html')

@app.route('/show_data')

def showData():
    # Retrieving uploaded file path from session
    file_list=os.listdir("./static/uploads")
    df = pd.DataFrame(columns=['file_name', 'S_sequence',"RBD_start","RBD_end","RBD_seq","Len","RBD_score_sum"])
    #fig
    wuhan_sum=1763.09
    wuhan_score=[8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161, 8.77161]
 
    position=list(range(331,532))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=position, y=wuhan_score, name = f'Wuhan:{wuhan_sum}',line=dict(color='royalblue', width=4)))
    for i in range(len(file_list)):
        if file_list[i].endswith(".fasta"):
            path= os.path.join("./static/uploads/",file_list[i])
            record = SeqIO.read(path, "fasta")
            record_aa=record.translate()
            res_start = re.search(r"NITNL", str(record_aa.seq)).start()
            res_end= re.search(r"GPKKST", str(record_aa.seq)).end()

            input_seq=record_aa[res_start:res_end]
            input_seq=str(input_seq.seq)
            df.loc[i]= file_list[i], str(record_aa.seq), res_start, res_end, input_seq, len(input_seq), ""
    for i in range(len(df)):
        input_seq=list(df["RBD_seq"][i])
        input_map=list([*map(col.get, input_seq)])
        input_score=[]
        for j in range(len(input_map)):
            bind=wuhan_m[j,input_map[j]]
            input_score.append(bind)
        input_sum=sum(input_score)
        df["RBD_score_sum"][i]= input_sum
        fig.add_trace(go.Scatter(x=position, y=input_score, name = f'{df["file_name"][i]}:{input_sum:.2f}'))
    #output
    df.to_csv("./static/outputs/summary_bywuhan.csv")
    fig.update_layout(title='RBD binding score',xaxis_title='Positions a.a.', yaxis_title='RBD score')
    fig.write_html("./static/outputs/score.html")
 
    # pandas dataframe to html table flask
    df.style.format(precision=2,formatter={('RBD_score_sum'): "{:.2f}"})
    df.style.set_properties(subset=['text'], **{'width': '200px'})
    uploaded_df_html = df.style.set_table_attributes('class="data"').to_html()
    # output plotly
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return render_template('result.html', data = uploaded_df_html, fig=graphJSON)
  
if __name__ == '__main__':
    from waitress import serve
    serve(app, host="0.0.0.0", port=8080)
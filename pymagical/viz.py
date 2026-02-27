import pandas as pd
import numpy as np
import re
import os
import json
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template

def parse_magical_output(file_path):
    """Parse the standard pymagical/MAGICAL tab-delimited output file."""
    data = []
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip('\n').split('\t')
            if len(parts) < 8:
                continue
            
            gene = parts[0]
            gene_chr = parts[1]
            gene_tss = parts[2]
            peak_chr = parts[3]
            peak_start = parts[4]
            peak_end = parts[5]
            peak_gene_prob = float(parts[6])
            tfs_raw = parts[7].strip().rstrip(',')
            
            # Extract TF details
            # Format: TF (prob, effect [L,B])
            # Example: ZNF557 (0.98, + [+,+])
            tf_matches = re.findall(r'([A-Za-z0-9_-]+)\s*\(([\d\.]+),\s*([\+\-])\s*\[([\+\-]),([\+\-])\]\)', tfs_raw)
            
            for tf_name, tf_prob, effect, l_dir, b_dir in tf_matches:
                data.append({
                    'Gene': gene,
                    'Gene_Chr': gene_chr,
                    'Gene_TSS': int(gene_tss),
                    'Peak_Chr': peak_chr,
                    'Peak_Start': int(peak_start),
                    'Peak_End': int(peak_end),
                    'Peak_Gene_Prob': peak_gene_prob,
                    'TF': tf_name,
                    'TF_Prob': float(tf_prob),
                    'Effect': effect,
                    'Looping_Dir': l_dir,
                    'Binding_Dir': b_dir,
                    'Circuit_Score': peak_gene_prob * float(tf_prob)
                })
    
    return pd.DataFrame(data)

def generate_report(input_file, output_html):
    """Generate a stylized HTML report from MAGICAL results."""
    df = parse_magical_output(input_file)
    
    if df.empty:
        print(f"No valid circuits found in {input_file}")
        return

    # --- Calculations ---
    total_circuits = len(df)
    unique_genes = df['Gene'].nunique()
    unique_tfs = df['TF'].nunique()
    
    # TF Frequency
    tf_counts = df['TF'].value_counts().reset_index()
    tf_counts.columns = ['TF', 'Count']
    
    # TF Effect Distribution
    tf_effects = df.groupby(['TF', 'Effect']).size().unstack(fill_value=0).reset_index()
    if '+' not in tf_effects.columns: tf_effects['+'] = 0
    if '-' not in tf_effects.columns: tf_effects['-'] = 0
    tf_effects['Total'] = tf_effects['+'] + tf_effects['-']
    tf_effects = tf_effects.sort_values('Total', ascending=False).head(20)

    # --- Plots ---
    # 1. TF Importance (Top 20)
    fig_tf = px.bar(tf_effects, x='TF', y=['+', '-'], 
                    title="Top 20 TFs by Circuit Involvement",
                    labels={'value': 'Number of Circuits', 'variable': 'Regulatory Effect'},
                    color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                    template="plotly_dark")
    fig_tf.update_layout(barmode='stack', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    
    # 2. Probability Distribution
    fig_dist = px.histogram(df, x='Circuit_Score', color='Effect',
                            title="Distribution of Combined Circuit Confidence",
                            labels={'Circuit_Score': 'Confidence (Peak-Gene Prob * TF-Peak Prob)'},
                            color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                            marginal="box", template="plotly_dark")
    fig_dist.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

    # 3. Peak-TSS Distance
    df['Distance'] = (df['Peak_Start'] + df['Peak_End'])/2 - df['Gene_TSS']
    df['Distance_kb'] = df['Distance'] / 1000
    fig_dist_tss = px.scatter(df, x='Distance_kb', y='TF_Prob', color='Effect',
                              hover_data=['Gene', 'TF'],
                              title="TF Confidence vs. Distance to TSS",
                              labels={'Distance_kb': 'Distance to TSS (kb)', 'TF_Prob': 'TF Binding Probability'},
                              color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                              template="plotly_dark")
    fig_dist_tss.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

    # Convert plots to JSON
    plots_json = {
        'tf_bar': fig_tf.to_json(),
        'score_dist': fig_dist.to_json(),
        'tss_dist': fig_dist_tss.to_json()
    }

    # --- HTML Template ---
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MAGICAL Results Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/dataTables.bootstrap5.min.css">
    <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    <style>
        :root {
            --bg-color: #0d1117;
            --card-bg: #161b22;
            --border-color: #30363d;
            --text-color: #c9d1d9;
            --accent-color: #58a6ff;
            --pos-color: #3fb950;
            --neg-color: #f85149;
        }
        body {
            background-color: var(--bg-color);
            color: var(--text-color);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
            margin: 0;
            padding: 20px;
        }
        .header {
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 20px;
            margin-bottom: 30px;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        .stat-card {
            background-color: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 20px;
            text-align: center;
        }
        .stat-value {
            font-size: 24px;
            font-weight: bold;
            color: var(--accent-color);
        }
        .plot-container {
            background-color: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 15px;
            margin-bottom: 20px;
        }
        .table-container {
            background-color: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 6px;
            padding: 20px;
            margin-top: 30px;
        }
        table.dataTable {
            background-color: var(--card-bg) !important;
            color: var(--text-color) !important;
        }
        .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, 
        .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing, 
        .dataTables_wrapper .dataTables_paginate {
            color: var(--text-color) !important;
        }
        .tag-pos { color: var(--pos-color); font-weight: bold; }
        .tag-neg { color: var(--neg-color); font-weight: bold; }
        input {
            background-color: var(--bg-color);
            border: 1px solid var(--border-color);
            color: var(--text-color);
            padding: 5px;
            border-radius: 4px;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>pymagical Inference Report</h1>
        <p>Input Source: <code>{{ input_file }}</code></p>
    </div>

    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-label">Total Circuits</div>
            <div class="stat-value">{{ total_circuits }}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Regulated Genes</div>
            <div class="stat-value">{{ unique_genes }}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Active TFs</div>
            <div class="stat-value">{{ unique_tfs }}</div>
        </div>
    </div>

    <div style="display: flex; flex-wrap: wrap; gap: 20px;">
        <div class="plot-container" style="flex: 1; min-width: 450px;" id="tf_bar"></div>
        <div class="plot-container" style="flex: 1; min-width: 450px;" id="score_dist"></div>
    </div>
    
    <div class="plot-container" id="tss_dist"></div>

    <div class="table-container">
        <h3>Circuit Detail Table</h3>
        <table id="resultsTable" class="display" style="width:100%">
            <thead>
                <tr>
                    <th>Gene</th>
                    <th>TF</th>
                    <th>Effect</th>
                    <th>Combined Prob</th>
                    <th>Peak Loc</th>
                    <th>Dist to TSS (kb)</th>
                    <th>L Dir</th>
                    <th>B Dir</th>
                </tr>
            </thead>
            <tbody>
                {% for row in table_data %}
                <tr>
                    <td>{{ row.Gene }}</td>
                    <td>{{ row.TF }}</td>
                    <td><span class="tag-{{ 'pos' if row.Effect == '+' else 'neg' }}">{{ row.Effect }}</span></td>
                    <td>{{ "%.4f"|format(row.Circuit_Score) }}</td>
                    <td>{{ row.Peak_Chr }}:{{ row.Peak_Start }}-{{ row.Peak_End }}</td>
                    <td>{{ "%.1f"|format(row.Distance_kb) }}</td>
                    <td>{{ row.Looping_Dir }}</td>
                    <td>{{ row.Binding_Dir }}</td>
                </tr>
                {% endfor %}
            </tbody>
        </table>
    </div>

    <script>
        const plots = {{ plots_json | safe }};
        Plotly.newPlot('tf_bar', JSON.parse(plots.tf_bar).data, JSON.parse(plots.tf_bar).layout);
        Plotly.newPlot('score_dist', JSON.parse(plots.score_dist).data, JSON.parse(plots.score_dist).layout);
        Plotly.newPlot('tss_dist', JSON.parse(plots.tss_dist).data, JSON.parse(plots.tss_dist).layout);

        $(document).ready(function() {
            $('#resultsTable').DataTable({
                "pageLength": 25,
                "order": [[ 3, "desc" ]],
                "dom": '<"top"f>rt<"bottom"lip><"clear">',
                "language": {
                    "search": "Filter results:"
                }
            });
        });
    </script>
</body>
</html>
    """

    # Prepare template data
    table_data = df.to_dict(orient='records')
    t = Template(html_template)
    rendered_html = t.render(
        input_file=os.path.basename(input_file),
        total_circuits=total_circuits,
        unique_genes=unique_genes,
        unique_tfs=unique_tfs,
        table_data=table_data,
        plots_json=json.dumps(plots_json)
    )

    with open(output_html, 'w') as f:
        f.write(rendered_html)
    
    print(f"Stylized report generated at: {output_html}")

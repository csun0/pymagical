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
    if not os.path.exists(file_path):
        return pd.DataFrame()
        
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
            # Supports both:
            # Python: TF (prob, effect [L,B])
            # MATLAB: TF (prob),
            tfs_list = [t.strip() for t in tfs_raw.split(',') if t.strip()]
            for tf_item in tfs_list:
                # 1. Try full Python format first
                full_match = re.match(r'([A-Za-z0-9_-]+)\s*\(([\d\.]+),\s*([\+\-])\s*\[([\+\-]),([\+\-])\]\)', tf_item)
                if full_match:
                    tf_name, tf_prob, effect, l_dir, b_dir = full_match.groups()
                else:
                    # 2. Try MATLAB format
                    ml_match = re.match(r'([A-Za-z0-9_-]+)\s*\(([\d\.]+)\)', tf_item)
                    if ml_match:
                        tf_name, tf_prob = ml_match.groups()
                        effect, l_dir, b_dir = '?', '?', '?'
                    else:
                        continue

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
    
    # TF Frequency for Bar Chart
    tf_effects = df.groupby(['TF', 'Effect']).size().unstack(fill_value=0).reset_index()
    if '+' not in tf_effects.columns: tf_effects['+'] = 0
    if '-' not in tf_effects.columns: tf_effects['-'] = 0
    tf_effects['Total'] = tf_effects['+'] + tf_effects['-']
    tf_effects_top = tf_effects.sort_values('Total', ascending=False).head(20)

    # TF Metrics for Bubble Plot
    # X: Unique Genes per TF, Y: Mean Circuit Score, Size: Total Circuits
    tf_bubble_data = df.groupby('TF').agg({
        'Gene': 'nunique',
        'Circuit_Score': 'mean',
        'Effect': lambda x: x.mode()[0] if not x.empty else '+'
    }).reset_index()
    tf_bubble_data.columns = ['TF', 'Unique_Genes', 'Mean_Confidence', 'Primary_Effect']
    # Add count for bubble size
    tf_counts_all = df['TF'].value_counts().reset_index()
    tf_counts_all.columns = ['TF', 'Circuit_Count']
    tf_bubble_data = tf_bubble_data.merge(tf_counts_all, on='TF')

    # --- Plots ---
    # 1. TF Importance (Top 20)
    fig_tf = px.bar(tf_effects_top, x='TF', y=['+', '-'], 
                    title="Top 20 TFs by Circuit Involvement",
                    labels={'value': 'Number of Circuits', 'variable': 'Effect'},
                    color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                    template="plotly_dark")
    fig_tf.update_layout(barmode='stack', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    
    # 2. Probability Distribution
    fig_dist = px.histogram(df, x='Circuit_Score', color='Effect',
                            title="Combined Circuit Confidence Distribution",
                            labels={'Circuit_Score': 'Combined Probability (P_peak_gene * P_tf_peak)'},
                            color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                            marginal="rug", template="plotly_dark")
    fig_dist.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

    # 3. Peak-TSS Distance
    df['Distance_kb'] = ((df['Peak_Start'] + df['Peak_End'])/2 - df['Gene_TSS']) / 1000
    fig_dist_tss = px.scatter(df, x='Distance_kb', y='Circuit_Score', color='Effect',
                              hover_data=['Gene', 'TF'],
                              title="Circuit Confidence vs. Distance to TSS",
                              labels={'Distance_kb': 'Distance to TSS (kb)', 'Circuit_Score': 'Combined Probability'},
                              color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                              template="plotly_dark")
    fig_dist_tss.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

    # 4. TF Bubble Plot (Influence Landscape)
    fig_bubble = px.scatter(tf_bubble_data, x='Unique_Genes', y='Mean_Confidence',
                            size='Circuit_Count', color='Primary_Effect',
                            hover_name='TF', text='TF',
                            title="TF Regulatory Influence Landscape",
                            labels={'Unique_Genes': 'Number of Unique Genes Regulated', 
                                    'Mean_Confidence': 'Mean Circuit Confidence'},
                            color_discrete_map={'+': '#3fb950', '-': '#f85149'},
                            size_max=40, template="plotly_dark")
    fig_bubble.update_traces(textposition='top center')
    fig_bubble.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')

    # Convert plots to JSON
    # Use engine='json' to ensure standard JSON instead of base64 binary blocks (bdata) 
    # which Plotly 6.0 uses but older JS versions may fail to parse.
    plots_json = {
        'tf_bar': fig_tf.to_json(engine='json'),
        'score_dist': fig_dist.to_json(engine='json'),
        'tss_dist': fig_dist_tss.to_json(engine='json'),
        'tf_bubble': fig_bubble.to_json(engine='json')
    }

    # --- HTML Template ---
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MAGICAL Results Report</title>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
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
            --text-muted: #8b949e;
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
            line-height: 1.5;
        }
        .header {
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 20px;
            margin-bottom: 30px;
        }
        .section-title {
            margin-top: 40px;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 1px solid var(--border-color);
            color: var(--accent-color);
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
            font-size: 28px;
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
        .doc-block {
            background-color: #1c2128;
            border-left: 4px solid var(--accent-color);
            padding: 15px;
            margin-bottom: 20px;
            font-size: 14px;
        }
        .takeaway {
            color: var(--text-muted);
            font-style: italic;
            margin-top: 10px;
            display: block;
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
        .tag-pos { color: var(--pos-color); font-weight: bold; }
        .tag-neg { color: var(--neg-color); font-weight: bold; }
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

    <h2 class="section-title">Global TF Analysis</h2>
    
    <div class="plot-container" id="tf_bubble"></div>
    <div class="doc-block">
        <strong>Plot: TF Regulatory Influence Landscape</strong><br>
        <strong>Calculation:</strong> X-axis counts unique gene symbols associated with each TF across all high-confidence circuits. Y-axis is the arithmetic mean of the 'Combined Probability' for all circuits involving that TF. Bubble size reflects total circuit count (one TF can regulate one gene via multiple peaks).<br>
        <span class="takeaway">Takeaway: TFs in the top-right are your "Master Regulators"—they control many genes with high statistical certainty. Large bubbles at the bottom-left may indicate TFs with broad but weak or noisy regulatory signatures.</span>
    </div>

    <div class="plot-container" id="tf_bar"></div>
    <div class="doc-block">
        <strong>Plot: TF Circuit Involvement</strong><br>
        <strong>Calculation:</strong> Raw count of TF occurrences in the output file, stacked by overall regulatory effect sign (+/-).<br>
        <span class="takeaway">Takeaway: Identifies the most frequently recovered TFs. A strong bias toward one color suggests a TF acting primarily as a global activator or repressor in this cell type.</span>
    </div>

    <h2 class="section-title">Circuit Mechanics & Confidence</h2>

    <div style="display: flex; flex-wrap: wrap; gap: 20px;">
        <div style="flex: 1; min-width: 450px;">
            <div class="plot-container" id="score_dist"></div>
            <div class="doc-block">
                <strong>Plot: Confidence Distribution</strong><br>
                <strong>Calculation:</strong> Combined Probability = (Peak-Gene Looping Prob) × (TF-Peak Binding Prob).<br>
                <span class="takeaway">Takeaway: A right-skewed distribution indicates a very clean, high-confidence signal. Multiple peaks might suggest distinct subpopulations of circuits with varying strengths.</span>
            </div>
        </div>
        <div style="flex: 1; min-width: 450px;">
            <div class="plot-container" id="tss_dist"></div>
            <div class="doc-block">
                <strong>Plot: Confidence vs. Distance</strong><br>
                <strong>Calculation:</strong> Distance = (Peak Midpoint) - (Gene TSS).<br>
                <span class="takeaway">Takeaway: Visualizes the spatial distribution of regulation. High-confidence points far from the TSS (e.g., >50kb) highlight important distal enhancer-looping events that standard methods often miss.</span>
            </div>
        </div>
    </div>

    <h2 class="section-title">Inferred Circuits Detail</h2>
    <div class="doc-block">
        <strong>Column Documentation:</strong>
        <ul>
            <li><strong>Combined Prob:</strong> The final circuit score ($P_{looping} \times P_{binding}$). Reflects the probability that the entire TF→Peak→Gene chain is active.</li>
            <li><strong>Dist to TSS:</strong> Distance in kilobases from the peak center to the gene's Transcription Start Site.</li>
            <li><strong>L/B Dir:</strong> The sign of the continuous regression weights for Looping (Peak-Gene) and Binding (TF-Peak).</li>
            <li><strong>Effect:</strong> The biological product of the weights. (+): Activator, (-): Repressor.</li>
        </ul>
    </div>

    <div class="table-container">
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
        Plotly.newPlot('tf_bubble', JSON.parse(plots.tf_bubble).data, JSON.parse(plots.tf_bubble).layout);

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

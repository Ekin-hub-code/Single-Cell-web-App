import zipfile
import os
import matplotlib.pyplot as plt
import scanpy as sc
from shiny import App, ui, render, reactive
from pathlib import Path
from sklearn.datasets import load_iris

extract_path = "extracted_files"

def extractzip(zip_path):
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
            print(f"Dosya {extract_path} dizinine çıkarıldı.")
    except Exception as e:
        return f"ZIP açılırken hata oluştu: {e}"

def load_h5ad(file_path):
    adata = sc.read(file_path)
    print("Veri başarıyla yüklendi.")
    return adata

app_ui = ui.page_fluid(
    ui.h2("Single cell viewer app"),
    ui.input_file("file_upload", "ZIP dosyası yükle", multiple=False, accept=[".zip"]),
    ui.output_plot("umap_plot"),
    ui.output_text_verbatim("status"),
    ui.input_action_button("action_button", "Process"),
)

def server(input, output, session):
    @reactive.event(input.action_button)
    def process_file():
        file_info = input.file_upload()
        if not file_info:
            return "No file uploaded."

        # ZIP dosyasını çıkar
        extractzip(file_info["datapath"])
        
        # .h5ad dosyasını bul ve yükle
        h5ad_files = list(Path(extract_path).glob("*.h5ad"))
        if not h5ad_files:
            return "No .h5ad file found."

        adata = load_h5ad(h5ad_files[0])

        # UMAP işlemleri
        sc.pp.pca(adata)
        sc.pp.neighbors(adata, use_rep="X_pca")
        sc.tl.umap(adata)
        sc.tl.leiden(adata)

        # Grafik çizimi
        plt.figure(figsize=(6, 4))
        sc.pl.umap(adata, color="Cell type", cmap="tab10", alpha=0.5, show=False)
        return plt.gcf()

    @output
    @render.text
    def status():
        return process_file()

    @output
    @render.plot
    def umap_plot():
        return process_file()

app = App(app_ui, server)



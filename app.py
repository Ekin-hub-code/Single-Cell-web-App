import zipfile
import os
import matplotlib.pyplot as plt
from shiny import App, ui, render, reactive
import umap
import numpy as np
from skylearn.datasets import load_iris

 #ZIP dosyasını açma
def extractzip(zip_path):   
 try:
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_path)
        print(f"Dosya {extract_path} dizinine çıkarıldı.")
 except Exception as e:
            return f"ZIP açılırken hata oluştu: {e}"

def load_h5ad(file_path):
    # .h5ad dosyasını yükleme
    adata = sc.read(file_path)
    print("Veri başarıyla yüklendi.")
    return adata
h5ad_files = list(extract_path.glob("*.h5ad"))
if not h5ad_files:
return ".h5ad_files can not be found."


# UI Tanımlama
app_ui = ui.page_fluid(
    ui.h2("Single cell viwer app"),
    ui.input_file("file_upload", "ZIP dosyası yükle", multiple=False, accept=[".zip"]),
    ui.output_plot("umap_plot"),
    ui.output_text_verbatim("status"),
    ui.input_action_button("action_button", "Action"),  
    ui.output_text("counter"),
)
app_ui = ui.page_sidebar(  
    ui.sidebar("Sidebar", position="right", bg="#f8f8f8"),  
    "Main content",
)  
def server(input, output, session):
    pass
    @output
    @render.text()
    @reactive.event(input.action_button)
    def counter():
        return f"{input.action_button()}"
    @reactive.calc
    def process_file():
        file_info = input.file_upload()
        if not file_info:
            return None


        sc.pp.pca(adata)
        sc.pp.neighbors(adata, use_rep="X_pca")
        sc.tl.umap(adata)

        sc.tl.leiden(adata)
        
        plt.figure(figsize=(6, 4))
        sc.pl.umap(adata, color="Cell type", cmap="tab10", alpha=0.5, show=False)
        
        return plt.gcf()



#ShinyApp'i çalıştır
app = App(app_ui, server)

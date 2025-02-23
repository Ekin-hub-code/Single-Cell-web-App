import zipfile
import os
import matplotlib.pyplot as plt
from shiny import App, ui, render, reactive
import umap
import numpy as np
from skylearn.datasets import load_iris

def unzip_file(zip_path, extract_path):
    #ZIP dosyasını açma
try:
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_path)
        print(f"Dosya {extract_path} dizinine çıkarıldı.")
    return extract_path

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
    ui.output_text_verbatim("status")
)
def server(input, output, session):
    @reactive.calc
    def process_file():
        file_info = input.file_upload()
        if not file_info:
            return None
data = load_iris()
X = data.data
Y = data.target

#UMAP modelini oluştur
umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2)
#Veriyi 2 boyuta indirgeme
X_umap = umap_model.fit_transform(X)
#Sonuçları görselleştirme
plt.scatter(X_umap[:, 0], X_umap[:,1], c=y, cmap='viridis', edgecolors='k', s=50)
plt.colorbar()
plt.title('UMAP - Iris veri Kümesi Boyut İndirgeme')
plt.show()


app = App(app_ui, server)

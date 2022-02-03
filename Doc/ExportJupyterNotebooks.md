# How to export Jupyter notebooks to HTML

Very useful to share your analysis with colleagues.


## Generate HTML

In a shell:

```shell
jupyter nbconvert --config $P3_REPO/Doc/nbconvert_config_html_exportOutput.py --template HtmlCollapsibleCells MyNotebook.ipynb
```

This generates a folder `htmlExport` with the html file and relative png figures inside.


## Publish with GitHub Pages


### Move to Pages folder

Move the folder `htmlExport` in `docs/JupyterNotebooks/`.
It is also recommended to rename it with the same name of html file, to finally have: `docs/JupyterNotebooks/MyNotebook/`


### Share the link with anybody

The public link to acces the notebook is then:

```
https://paulscherrerinstitute.github.io/PSIPositronProduction/JupyterNotebooks/MyNotebook/MyNotebook.html
```

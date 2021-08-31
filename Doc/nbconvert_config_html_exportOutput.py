c = get_config()

# print(c.NbConvertApp.notebooks)
# print(c.NbConvertApp.output_base)
c.NbConvertApp.export_format = 'html'
# Output directory
c.FilesWriter.build_directory = 'htmlExport'
c.NbConvertApp.output_files_dir = ''
# Filename of images
# keep it short, too long paths do not work in GitLab Pages
# Default filename
# ExtractOutputPreprocessor.output_filename_template = '{unique_key}_{cell_index}_{index}{extension}'
# Custom filename
c.ExtractOutputPreprocessor.output_filename_template = 'Out_{cell_index}_{index}{extension}'
# Add additional preprocessors
c.HTMLExporter.preprocessors = [
    'nbconvert.preprocessors.ExtractOutputPreprocessor',
]
# Default preprocessors
# c.HTMLExporter.default_preprocessors = [
#     'nbconvert.preprocessors.TagRemovePreprocessor',
#     'nbconvert.preprocessors.RegexRemovePreprocessor',
#     'nbconvert.preprocessors.coalesce_streams',
#     'nbconvert.preprocessors.CSSHTMLHeaderPreprocessor',
#     'nbconvert.preprocessors.HighlightMagicsPreprocessor',
# ]

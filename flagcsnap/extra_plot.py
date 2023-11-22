import pandas as pd

class ExtraPlotter:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):       

        results = self.results
        den_data = self.den_data
        dendrogram = self.dendrogram
        results_ncfeatures = self.results_ncfeatures

        results["temp"]= results["assembly_accession"].str.split(".").str[0]
        results=pd.merge(results,self.type_strains[["bacdive_id","strain_number_header","type_strain_assembly"]],left_on="temp",right_on="type_strain_assembly",how="left").drop(columns=["temp","type_strain_assembly"])
        den_data["temp"]= den_data["assembly_accession"].str.split(".").str[0]
        den_data = pd.merge(den_data,self.type_strains[["bacdive_id","strain_number_header","type_strain_assembly"]],left_on="temp",right_on="type_strain_assembly",how="left").drop(columns=["temp"])
        den_data['bckg_color'] = den_data['strain_number_header'].apply(
        lambda x: [255, 255, 255, 0] if pd.isnull(x) else [255, 228, 184, 255])    
        data=results.to_json(orient="records")
        dend_data=den_data.to_json(orient="records")
        dend=dendrogram.to_json(orient="records")
        nc_data=results_ncfeatures.to_json(orient="records")


        first_part=""" 
        <!DOCTYPE html>
        <html>
        <head>
        <title>Deck.gl Orthographic View</title>
        <script src="https://unpkg.com/deck.gl@^8.0.0/dist.min.js"></script>
        <script>
            data={data}
            dendrogram={dend}
            dend_data={dend_data}
            nc_data={nc_data}    
        </script>
        </head>
        """
        first_part=first_part.format(data=data,dend=dend,dend_data=dend_data,nc_data=nc_data)
        second_part=""" 
        <body>
        <div id="container" style="width:100%; height:100%; position:absolute;"></div>

        <script type="text/javascript">
            // // Define URL
            // const DATA_URL = './data.json';

            // // Fetch data from URL
            // fetch(DATA_URL)
            //   .then(response => response.json())
            //   .then(data => {
                // Define view and layers
                const {DeckGL, PolygonLayer, OrthographicView, TextLayer, PathLayer, LineLayer} = deck;

                const deckgl = new DeckGL({
                container: 'container',
                width: '100%',
                height: '100%',
                views: [new OrthographicView({id: 'ortho', controller: true})],
                initialViewState: {
                    target: [0, 0, 0],
                    zoom: 1
                },
                layers: [
                new LineLayer({
                    id: 'line-layer',
                    data: dend_data, // Assuming that this is where your data is
                    getSourcePosition: d => d.start_line,
                    getTargetPosition: d => d.end_line,
                    getWidth: 1,
                    getColor: [125,125,125,255]
                }),
                    new PolygonLayer({
                    id: 'polygon-layer',
                    data: data,
                    pickable: true,
                    autoHighlight: true,
                    stroked: true,
                    filled: true,
                    wireframe: true,
                    lineWidthMinPixels: 1,
                    getPolygon: d => d.coordinates,
                    getFillColor: d => d.fillcolor,  // Use the fillcolor property from your data
                    getLineColor: d => d.linecolor,  // Black color
                    getLineWidth: d => d.linetype, // Set line width
                    onClick: ({object}) => {
                        if (object && object.id) {
                            const originalProteinId = object.id;
                            let proteinId;
                            let link;
                    
                            if (originalProteinId.startsWith('IMGVR_')) {
                                proteinId = originalProteinId.split('|')[0];
                                link = `https://img.jgi.doe.gov/cgi-bin/vr/main.cgi?section=ViralBrowse&page=uviginfo&uvig_id=${proteinId}`;
                            } else if (originalProteinId.startsWith('IMGPR_')) {
                                proteinId = originalProteinId.split('|')[0];
                                link = `https://img.jgi.doe.gov/cgi-bin/plasmid/main.cgi?section=PlasmidBrowse&page=plasmidinfo&plasmid_id=IMGPR_plasmid_3300053491_061425`;
                            } else {
                                proteinId = originalProteinId;
                                link = `https://www.ncbi.nlm.nih.gov/protein/${proteinId}`;
                            }
                    
                            window.open(link, '_blank');
                        }
                    }  
                    }),
                    new PathLayer({
                    id: 'path-layer',
                    data: dendrogram, // Assuming that this is where your data is
                    pickable: false,
                    widthScale: 0.5,
                    widthMinPixels: 0,
                    getPath: d => d.path,
                    getWidth: 1,
                    opacity: 1,
                    getColor: [125,125,125,255]
                    }),
                    new TextLayer({
                    id: 'text-layer',
                    data: dend_data,  // assuming that this is where your data is
                    pickable: false,
                    getPosition: d => d.coordinates,
                    getText: d => d.species,
                    getColor: [0, 0, 0],
                    getAngle: 0,
                    getSize: 7.5,
                    background: true,
                    getBackgroundColor: d => d.bckg_color, 
                    getTextAnchor: 'start',
                    getAlignmentBaseline: 'center',
                    sizeUnits: 'meters',
                    fontFamily: 'Roboto, sans-serif',
                    sizeScale: 1,
                    }),
                    new PolygonLayer({
                    id: 'nc-layer',
                    data: nc_data,
                    pickable: true,
                    autoHighlight: true,
                    stroked: false,
                    filled: true,
                    wireframe: true,
                    lineWidthMinPixels: 1,
                    getPolygon: d => d.coordinates,
                    getFillColor: d => d.fillcolor,  // Use the fillcolor property from your data
                    onClick: ({object}) => {
                        if (object && object.id) {
                        const proteinId = object.id;
                        window.open(`https://www.ncbi.nlm.nih.gov/protein/${proteinId}`, '_blank');
                        }
                    }
                    })
                ],
                getTooltip: ({object}) => object && {
                    html: `
                    <div><b>SeqId:</b> ${object.seqid}</div>
                    <div><b>Start:</b> ${object.rel_start}</div>
                    <div><b>End:</b> ${object.rel_end}</div>
                    <div><b>Protein ID:</b> ${object.id}</div>
                    <div><b>Product:</b> ${object.product}</div>
                    <div><b>Padloc System:</b> ${object.padloc_system || 'N/A'}</div>
                    <div><b>Padloc Gene:</b> ${object.padloc_gene || 'N/A'}</div>
                    <div><b>Deffinder System:</b> ${object.deffinder_system || 'N/A'}</div>
                    <div><b>Deffinder Gene:</b> ${object.deffinder_gene || 'N/A'}</div>
                    <div><b>CCTyper System:</b> ${object.cctyper_system || 'N/A'}</div>
                    <div><b>CCTyper Gene:</b> ${object.cctyper_gene || 'N/A'}</div>
                    <div><b>BacDive ID:</b> ${object.bacdive_id || 'N/A'}</div>
                    <div><b>Type strain(s):</b> ${object.strain_number_header || 'N/A'}</div>
                    <div><b>Family ID:</b> <span style="color: rgba(${object.fillcolor}); font-weight: bold;">${object.fam_cluster || 'N/A'}</span></div>
                    <div><b>Non-coding feature:</b> <span style="color: rgba(${object.fillcolor}); font-weight: bold;">${object.nc_feature || 'N/A'}</span></div>
                    `,
                    style: {
                    fontFamily: 'Helvetica Neue, sans-serif',
                    color: 'black',
                    backgroundColor: '#fff',
                    padding: '5px',
                    border: '1px solid #ccc'
                    }
                }
                });
        </script>
        </body>
        </html>
        """
        f = open(self.output+"/"+self.output+"_extra.html","w+")
        f.write(first_part+second_part)
        f.close()

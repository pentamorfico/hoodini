import pandas as pd
import zlib,base64

class ExtraPlotter:
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run(self):       

        results = self.results
        den_data = self.den_data
        dendrogram = self.dendrogram
        if self.domains:
            domains = self.domains_data
        ticks = self.ticks
        tick_text = self.tick_text

        results["temp"]= results["assembly_accession"].str.split(".").str[0]
        results=pd.merge(results,self.type_strains[["bacdive_id","strain_number_header","type_strain_assembly"]],left_on="temp",right_on="type_strain_assembly",how="left").drop(columns=["temp","type_strain_assembly"])
        columns = ["seqid","id","strand","species","rel_start","rel_end","product","padloc_system","padloc_gene","deffinder_system","deffinder_gene","cctyper_system","cctyper_gene","bacdive_id","strain_number_header","fam_cluster","coordinates","fillcolor","linecolor","linetype","text_coordinates","prevalence","sequence","assembly_accession"]
        results = results[columns]
        den_data["temp"]= den_data["assembly_accession"].str.split(".").str[0]
        den_data = pd.merge(den_data,self.type_strains[["bacdive_id","strain_number_header","type_strain_assembly"]],left_on="temp",right_on="type_strain_assembly",how="left").drop(columns=["temp"])
        den_data['bckg_color'] = den_data['strain_number_header'].apply(
        lambda x: [255, 255, 255, 0] if pd.isnull(x) else [255, 228, 184, 255])   

        data = results.to_json(orient="records")
        data = zlib.compress(data.encode())
        data = base64.b64encode(data)
        data = data.decode()

        dend_data = den_data.to_json(orient="records")
        dend_data = zlib.compress(dend_data.encode())
        dend_data = base64.b64encode(dend_data)
        dend_data = dend_data.decode()

        dendrogram = dendrogram.to_json(orient="records")
        dendrogram = zlib.compress(dendrogram.encode())
        dendrogram = base64.b64encode(dendrogram)
        dendrogram = dendrogram.decode()

        if self.domains:
            domains = domains.to_json(orient="records")
            domains = zlib.compress(domains.encode())
            domains = base64.b64encode(domains)
            domains = domains.decode()
        else:
            domains = {}

        ticks = ticks.to_json(orient="records")
        ticks = zlib.compress(ticks.encode())
        ticks = base64.b64encode(ticks)
        ticks = ticks.decode()

        tick_text = tick_text.to_json(orient="records")
        tick_text = zlib.compress(tick_text.encode())
        tick_text = base64.b64encode(tick_text)
        tick_text = tick_text.decode()

        first_part=""" 
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title> Flagcsnap</title>
            <script src="https://unpkg.com/react@17/umd/react.development.js"></script>
            <script src="https://unpkg.com/react-dom@17/umd/react-dom.development.js"></script>
            <script src="https://unpkg.com/pako@2.1.0/dist/pako.min.js"></script>
            

            <!-- Babel for JSX -->
            <script src="https://unpkg.com/@babel/standalone@7/babel.min.js"></script>

            <!-- Deck.gl and Mapbox -->
            <script src="https://unpkg.com/deck.gl@8/dist.min.js"></script>

            <script src="https://unpkg.com/@mui/material@latest/umd/material-ui.development.js"></script>

            <script src="https://cdn.rawgit.com/arose/ngl/v2.0.0-dev.31/dist/ngl.js"></script>
            <!-- Fonts -->
            <link rel="preconnect" href="https://fonts.googleapis.com">
            <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
            <link href="https://fonts.googleapis.com/css2?family=Roboto&display=swap" rel="stylesheet">
            <link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto:300,400,500,700&display=swap" />
            <link rel="stylesheet" href="https://fonts.googleapis.com/icon?family=Material+Icons" />

            <script type="text/javascript" src="https://threejs.org/build/three.js"></script>

            <link href="https://fonts.cdnfonts.com/css/helvetica-neue-5" rel="stylesheet">

            <script>

            function decompressData(compressedBase64) {{
                var compressedBytes = atob(compressedBase64);
                var charData = compressedBytes.split('').map(function(x) {{ return x.charCodeAt(0); }});
                var binData = new Uint8Array(charData);
                return pako.inflate(binData, {{ to: 'string' }});
            }}

            data="{data}";
            dendrogram="{dendrogram}";
            dend_data="{dend_data}";
            ticks="{ticks}";
            tick_text="{tick_text}";
            domains="{domains}";


            data = JSON.parse(decompressData(data));
            dendrogram = JSON.parse(decompressData(dendrogram));
            dend_data = JSON.parse(decompressData(dend_data));
            ticks = JSON.parse(decompressData(ticks));
            tick_text = JSON.parse(decompressData(tick_text));
            domains = JSON.parse(decompressData(domains));

            </script>
            <style>
                .info-panel {{
                    background: rgba(255, 255, 255, 0.312);
                    border-radius: 2px;
                    box-shadow: 0 4px 30px rgba(0, 0, 0, 0.1);
                    backdrop-filter: blur(2px);
                    -webkit-backdrop-filter: blur(2px);
                    overflow: auto;
                    word-wrap: break-word;
                    position: absolute;
                    cursor: grab;
                    top: 50px;
                    left: 50px;
                    max-width: 300px;
                    padding: 30px;
                    border-radius: 5px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.3);
                    z-index: 9;
                    font-family: 'Roboto', sans-serif;
                    color: rgba(0, 0, 0, 1);
                    overflow-y: auto;
                    max-height: 90%;
                }}

                .info-panel a {{
                    color: #0000EE;
                    text-decoration: none;
                }}

                .info-panel a:hover {{
                    text-decoration: underline;
                }}

                .info-panel p {{
                    margin: 5px 0;
                }}

                .info-panel strong {{
                    font-weight: bold;
                }}
            </style>
        </head>
        """

        first_part=first_part.format(data=data,dendrogram=dendrogram,dend_data=dend_data,domains=domains,ticks=ticks,tick_text=tick_text)
        
        second_part=""" 
        <body>
        <div id="root" style="width:100%; height:100%; position:absolute;"></div>

        <script type="text/babel">
                const { useEffect, useRef, useState } = React;
                const { DeckGL, PolygonLayer, OrthographicView, LineLayer, TextLayer, PathLayer, OrthographicController } = deck;

                function App() {
                    const [selectedObject, setSelectedObject] = useState(null);
                    const [key, setKey] = useState(0); // Inicializa un contador de clave

                    const handleObjectSelect = (object) => {
                        setSelectedObject(object);
                        setKey(prevKey => prevKey + 1); // Incrementa la clave
                    };

                    const handleClose = () => {
                        setSelectedObject(null); // Esto ocultará el SidePanel
                    };




                    return (
                        <div style={{ width: '100%', height: '100%' }}>
                            <DeckGLReactComponent onObjectSelect={handleObjectSelect} selectedObject={selectedObject} />
                            {selectedObject && <SidePanel key={key} object={selectedObject} onClose={handleClose} />}
                        </div>
                    );
                }


                function DeckGLReactComponent({ onObjectSelect }) {
                    const [selectedPolygonId, setSelectedPolygonId] = useState(null);
                    const containerRef = useRef(null);

                    useEffect(() => {



                        class MyOrthographicController extends OrthographicController {


                            handleKeyDown(event) {
                                if (event.key === 'Control') {
                                    this.ctrlKeyPressed = true;
                                }
                            }

                            handleKeyUp(event) {
                                if (event.key === 'Control') {
                                    this.ctrlKeyPressed = false;
                                }
                            }


                            _onWheel(event) {
                                if (!this.scrollZoom && !this.dragPan) {
                                    return false;
                                }

                                const controllerState = this.controllerState;
                                let newControllerState;
                                const interactionState = {};

                                // Check if the Control key is pressed
                                const isControlPressed = event.srcEvent.ctrlKey;
                                const isAltPressed = event.srcEvent.AltKey;


                                if (isControlPressed && !this.ctrlKeyPressed){


                                    const pos = this.getCenter(event);
                                    if (!this.isPointInBounds(pos, event)) {
                                        return false;
                                    }
                                    event.srcEvent.preventDefault();

                                    const { speed = 0.05, smooth = false } = this.scrollZoom === true ? {} : this.scrollZoom;
                                    const { delta } = event;

                                    // Map wheel delta to relative scale
                                    let scale = 2 / (1 + Math.exp(-Math.abs(delta * speed)));
                                    if (delta < 0 && scale !== 0) {
                                        scale = 1 / scale;
                                    }

                                    const newControllerState = this.controllerState.zoom({ pos, scale });
                                    this.updateViewport(
                                        newControllerState,
                                        { ...this._getTransitionProps({ around: pos }), transitionDuration: smooth ? 250 : 1 },
                                        {
                                            isZooming: true,
                                            isPanning: true
                                        }
                                    );
                                    return true;

                                }
                                
                                else if (isControlPressed) {

                                    const pos = this.getCenter(event);
                                    if (!this.isPointInBounds(pos, event)) {
                                        return false;
                                    }
                                    event.srcEvent.preventDefault();


                                    const { speed = 0.01, smooth = false } = this.scrollZoom === true ? {} : this.scrollZoom;
                                    const { delta } = event;

                                    // Map wheel delta to relative scale
                                    let scale = 2 / (1 + Math.exp(-Math.abs(delta * speed)));
                                    if (delta < 0 && scale !== 0) {
                                        scale = 1 / scale;
                                    }


                                    const newControllerState = this.controllerState.zoom({ pos, scale });
                                    this.updateViewport(
                                        newControllerState,
                                        { ...this._getTransitionProps({ around: pos }), transitionDuration: smooth ? 250 : 1 },
                                        {
                                            isZooming: true,
                                            isPanning: true
                                        }
                                    );
                                    return true;
                                } else {
                                const pos = this.getCenter(event);
                                if (!this.isPointInBounds(pos, event)) {
                                    return false;
                                }
                                event.srcEvent.preventDefault();

                                const deltaX = event.srcEvent.deltaX;
                                const deltaY = event.srcEvent.deltaY;
                                const moveSpeed = 150; // Adjust this value as needed
                                const transitionDuration = 200; // Duration of the transition in milliseconds
                                const easingFunction = t => t * (2 - t); // Simple ease-out function

                                let newControllerState = this.controllerState;

                                if (Math.abs(deltaX) > Math.abs(deltaY)) {
                                    // Horizontal movement only
                                    if (deltaX > 0) {
                                        newControllerState = newControllerState.moveLeft(moveSpeed);
                                    } else if (deltaX < 0) {
                                        newControllerState = newControllerState.moveRight(moveSpeed);
                                    }
                                } else {
                                    // Vertical movement only
                                    if (deltaY < 0) {
                                        newControllerState = newControllerState.moveDown(moveSpeed);
                                    } else if (deltaY > 0) {
                                        newControllerState = newControllerState.moveUp(moveSpeed);
                                    }
                                }

                                this.updateViewport(newControllerState, {
                                    transitionDuration: transitionDuration,
                                    transitionEasing: easingFunction
                                }, {
                                    isPanning: true
                                });

                                return true;

                            }

                                // If you have the original implementation of _onWheel, it should go here
                                // For now, let's return false as a placeholder
                                return false;
                            }
                            
                        }


                        let zoomMode = "all"
                        const deckgl = new DeckGL({
                            container: containerRef.current,
                            width: '100%',
                            height: '100%',
                            inertia: true,
                            views: [new OrthographicView({
                                id: "main",
                                controller: {
                                    type:MyOrthographicController, inertia: true,zoomAxis:"all"
                                },
                            }
                            )],
                            onViewStateChange: ({viewState}) => {
      
                                document.addEventListener('keydown', (event) => {
                                    if (event.key === 'Alt' || event.altKey) {
                                        event.preventDefault();

                                        // When the Alt key is pressed
                                        const nextZoomMode = 'X'; // Set to 'X' or your desired mode
                                        if (zoomMode !== nextZoomMode) {
                                            zoomMode = nextZoomMode;
                                            deckgl.setProps({
                                                controller: {type:MyOrthographicController, inertia: true, zoomAxis: "X" }
                                            });
                                        }
                                    }
                                });

                                document.addEventListener('keyup', (event) => {
                                    if (!event.altKey) {

                                        // When the Alt key is released
                                        const nextZoomMode = 'all'; // Set back to 'all' or your normal mode
                                        if (zoomMode !== nextZoomMode) {
                                            zoomMode = nextZoomMode;
                                            deckgl.setProps({
                                                controller: {type:MyOrthographicController, inertia: true, zoomAxis: zoomMode}
                                            });
                                        }
                                    }
                                });   

                                document.addEventListener('keydown', (event) => {
                                    if (event.shiftKey) {
                                        event.preventDefault(); // Prevent default action for Shift key
                                        const nextZoomMode = 'Y';
                                        if (zoomMode !== nextZoomMode) {
                                            zoomMode = nextZoomMode;
                                            deckgl.setProps({
                                                controller: {type:MyOrthographicController, inertia: true, zoomAxis: zoomMode}
                                            });
                                        }
                                    }
                                });

                                document.addEventListener('keyup', (event) => {
                                    if (!event.shiftKey && zoomMode === 'Y') {
                                        event.preventDefault(); // Prevent default action when releasing Shift key
                                        const nextZoomMode = 'all';
                                        if (zoomMode !== nextZoomMode) {
                                            zoomMode = nextZoomMode;
                                            deckgl.setProps({
                                                controller: {type:MyOrthographicController, inertia: true, zoomAxis: zoomMode}
                                            });
                                        }
                                    }
                                });

                                
                            },
                            initialViewState: {
                                target: [0, 0],
                                zoom: [0,0]
                            },
                        layers: [
                            new LineLayer({
                                id: 'track_lines',
                                data: dend_data,
                                getSourcePosition: d => d.start_line,
                                getTargetPosition: d => d.end_line,
                                getWidth: 1,
                                getColor: [125, 125, 125, 255]
                            }),
                            new LineLayer({
                                id: 'tick_lines',
                                data: ticks,
                                getSourcePosition: d => d.start,
                                getTargetPosition: d => d.end,
                                getWidth: 1,
                                getColor: [125, 125, 125, 255]
                            }),
                            new TextLayer({
                                id: 'tick_text_layer',
                                data: tick_text,  // assuming that this is where your data is
                                pickable: false,
                                getPosition: d => d.coordinates,
                                getText: d => d.text,
                                getColor: [0, 0, 0],
                                getAngle: 0,
                                getSize: 2.5,
                                background: false,
                                getTextAnchor: 'middle',
                                getAlignmentBaseline: 'top',
                                sizeUnits: 'meters',
                                fontFamily: 'Roboto, sans-serif',
                                sizeScale: 1,
                            }),
                            new PolygonLayer({
                                id: 'gene_layer',
                                data: data,
                                pickable: true,
                                autoHighlight: true,
                                stroked: true,
                                filled: true,
                                wireframe: true,
                                lineWidthMinPixels: 1,
                                getPolygon: d => d.coordinates,
                                getFillColor: d => d.fillcolor,  // Use the fillcolor property from your data
                                getLineWidth: d => d.linetype, // Set line width
                                getLineColor: d => d.linecolor,
                                onClick: (info, event) => {
                                    console.log('Clicked:', info, event);
                                    onObjectSelect(info.object);
                                }
                            }),
                            new PolygonLayer({
                                id: 'domain_layer',
                                data: domains,
                                pickable: true,
                                autoHighlight: true,
                                stroked: true,
                                filled: true,
                                wireframe: true,
                                lineWidthMinPixels: 1,
                                getPolygon: d => d.coordinates,
                                getFillColor: d => d.colors,
                                getLineWidth: 0.05, // Set line width
                            }),
                            new TextLayer({
                                id: 'genetext_layer',
                                data: data,  // assuming that this is where your data is
                                pickable: false,
                                getPosition: d => d.text_coordinates,
                                getText: d => d.fam_cluster,
                                getColor: [0, 0, 0],
                                getAngle: 0,
                                getSize: 2.5,
                                getColor: d => d.linecolor,
                                background: false,
                                getTextAnchor: 'middle',
                                getAlignmentBaseline: 'top',
                                sizeUnits: 'meters',
                                fontFamily: 'Roboto, sans-serif',
                                sizeScale: 1,
                            }),
                            new PathLayer({
                                id: 'dendrogram_layer',
                                data: dendrogram,
                                pickable: false,
                                widthScale: 0.5,
                                widthMinPixels: 0,
                                getPath: d => d.path,
                                getWidth: 1,
                                opacity: 1,
                                getColor: [125, 125, 125, 255]
                            }),
                            new TextLayer({
                                id: 'species_layer',
                                data: dend_data,
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
                        ],
                        getTooltip: ({object, layer}) => {
                            if (object && layer.id === "gene_layer") {
                                return {
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
                                    <div><b>Prevalence:</b> ${(object.prevalence * 100).toFixed(2)}%</div>
                                    `,
                                    style: {
                                        fontFamily: 'Roboto, sans-serif',
                                        color: 'black',
                                        backgroundColor: '#fff',
                                        padding: '5px',
                                        border: '1px solid #ccc'
                                    }
                                };
                            } else if (object && layer.id === "domain_layer") {
                                return {
                                    html: `
                                    <div><b>Domain:</b> ${object.domain_id}</div>
                                    <div><b>Short name:</b> ${object.short_name}</div>
                                    <div><b>Description:</b> ${object.description}</div>
                                    <div><b>Clan:</b> ${object.clan  || ''}</div>
                                    <div><b>E-value:</b> ${object.e_value}</div>
                                    <div><b>Bitscore:</b> ${object.bit_score || 'N/A'}</div>
                                    `,
                                    style: {
                                        fontFamily: 'Roboto, sans-serif',
                                        color: 'black',
                                        backgroundColor: '#fff',
                                        padding: '5px',
                                        border: '1px solid #ccc'
                                    }
                                };
                            }
                        }
                    });

                    return () => deckgl && deckgl.finalize();
                }, []);

                return (
                    <div ref={containerRef} style={{ width: '100%', height: '100%' }} />
                );
            }

            function SidePanel({ object, onClose }) {

            const [showContextSliders, setShowContextSliders] = useState(false);
            const [upstreamValue, setUpstreamValue] = useState(5000);
            const [downstreamValue, setDownstreamValue] = useState(5000);



            const [showWarning, setShowWarning] = useState(false);
            const [pdbData, setPdbData] = useState(null); // State to hold pdbData
            const [showStructure, setShowStructure] = useState(false);

            const { Divider, Typography,Box, Grid, Chip, IconButton, Paper, Fade, Button, Slider, Input } = MaterialUI;

            const [isDragging, setIsDragging] = useState(false);
            const [showProteinSequence, setShowProteinSequence] = useState(false);
            const [showTaxonomy, setShowTaxonomy] = useState(false);
            const [position, setPosition] = useState({ x: 50, y: 50 });
            const panelRef = useRef(null);

            const startDrag = (e) => {
            if (showStructure || showContextSliders ) return; // No permitir arrastrar si showStructure es true
            setIsDragging(true);
        };

            const onDrag = (e) => {
                if (!isDragging) return;
                const newPosition = {
                    x: position.x + e.movementX,
                    y: position.y + e.movementY,
                };
                setPosition(newPosition);
            };

            const stopDrag = () => {
                setIsDragging(false);
            };

            useEffect(() => {
                const handleMouseMove = (e) => {
                    onDrag(e);
                };

                const handleMouseUp = () => {
                    stopDrag();
                };

                if (isDragging) {
                    document.addEventListener('mousemove', handleMouseMove);
                    document.addEventListener('mouseup', handleMouseUp);
                }

                return () => {
                    document.removeEventListener('mousemove', handleMouseMove);
                    document.removeEventListener('mouseup', handleMouseUp);
                };
            }, [isDragging, onDrag]);
            
            
            async function fetchProteinStructure(sequence, refseq_id) {
                
            const esm_url = 'https://api.esmatlas.com/foldSequence/v1/pdb/';
            const uniprot_url = "https://rest.uniprot.org/uniprotkb/search?query=" + refseq_id;
            const afdb_url = "https://alphafold.ebi.ac.uk/api/prediction/";

            try {
                const uniprot_response = await fetch(uniprot_url);
                console.log(uniprot_response);
                if (!uniprot_response.ok) {
                    throw new Error(`HTTP error! status: ${uniprot_response.status}`);
                }

                const uniprot_data = await uniprot_response.json();
                
                if (uniprot_data.results.length === 0){
                    if (sequence.length > 400) {
                        setShowWarning(true);
                        
                        return;
                    } 

                    else {
                        setShowWarning(false);

                    try {
                        const response = await fetch(esm_url, {
                            method: 'POST',
                            headers: {
                                'Content-Type': 'text/plain'
                            },
                            body: sequence
                        });

                        if (!response.ok) {
                            throw new Error(`HTTP error! status: ${response.status}`);
                        }

                        const data = await response.text();
                        return data;
                    } catch (error) {
                        console.error('Error fetching protein structure:', error);
                    }
                }

                }

                try {
                    const results = uniprot_data.results;
                    for (const result of results) {
                        const uniProtKBCrossReferences = result.uniProtKBCrossReferences;
                        for (const db_entry of uniProtKBCrossReferences) {
                            if (db_entry.database === 'AlphaFoldDB') {
                                const af_id = db_entry.id;
                                console.log("AlphaFoldDB ID:", af_id);
                                const afdb_data = await (await fetch(afdb_url + af_id)).json();
                                const pdb_url = afdb_data[0].pdbUrl;
                                console.log(pdb_url);
                                const afdb_pdb = await (await fetch(pdb_url)).text();
                                return afdb_pdb;
                            }
                        }
                    }
                } catch (error) {
                    console.error('No AFDB-ID found:', error);
                }

            } catch (error) {
                console.error('Error fetching protein structure:', error);
            }
        }
            
            function handleShowStructure() {

            if (showStructure) {
                setShowStructure(false);
                return;
            }
            
            const proteinSequence = object.sequence;
            const refseq_id = object.id;
            fetchProteinStructure(proteinSequence, refseq_id).then(pdbData => {
                if (pdbData) {
                    setShowStructure(true);
                    console.log('pdbData:', pdbData);
                    setPdbData(pdbData);
                    // Aquí puedes usar pdbData para cargar en el visor NGL
                    var stage = new NGL.Stage("viewport" , {backgroundColor: "white"});
                    stage.loadFile(new Blob([pdbData], {type: 'text/plain'}), { ext: 'pdb'}).then(function(o) {
                o.addRepresentation("cartoon");
                o.autoView();
            }
            );        }
            });

        }

        const handlePdbDataDownload = () => {

                // if !pdbData show consolelog saying that no pdbdata

                if (!pdbData) {
                    console.log('No pdbData');
                    return;
                }

                // Create a Blob from the pdbData
                const blob = new Blob([pdbData], { type: 'text/plain' });
                const url = URL.createObjectURL(blob);

                // Create a link and trigger the download
                const a = document.createElement('a');
                a.href = url;
                //use the protein id to name the file
                a.download = `${object.id}.pdb`;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);

                // Revoke the URL
                URL.revokeObjectURL(url);
            };

            // Function to handle changes in sliders
            const handleSliderChange = (type, newValue) => {
                if (type === 'upstream') {
                    setUpstreamValue(newValue);
                } else {
                    setDownstreamValue(newValue);
                }
            };

            // Function to handle changes in input fields
            const handleInputChange = (event, type) => {
                const value = event.target.value === '' ? 0 : Number(event.target.value);
                if (type === 'upstream') {
                    setUpstreamValue(value > 20000 ? 20000 : value);
                } else {
                    setDownstreamValue(value > 20000 ? 20000 : value);
                }
            };

            // Function to handle blur event on input fields
            const handleBlur = (type) => {
                if (type === 'upstream' && (upstreamValue < 0 || upstreamValue > 20000)) {
                    setUpstreamValue(upstreamValue < 0 ? 0 : 20000);
                } else if (type === 'downstream' && (downstreamValue < 0 || downstreamValue > 20000)) {
                    setDownstreamValue(downstreamValue < 0 ? 0 : 20000);
                }
            };

            const handleDownload = async () => {
                let seqStart, seqEnd;

                if (object.flipped === "flipped") {
                    seqStart = object.end + upstreamValue;
                    seqEnd = object.start - downstreamValue;
                } else {
                    seqStart = object.start - upstreamValue;
                    seqEnd = object.end + downstreamValue;
                }

                const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${object.seqid}&rettype=gb&retmode=text&seq_start=${seqStart}&seq_stop=${seqEnd}`;

                try {
                    const response = await fetch(url);
                    if (!response.ok) {
                        throw new Error(`HTTP error! status: ${response.status}`);
                    }
                    const data = await response.text();

                    // Create a Blob from the data
                    const blob = new Blob([data], { type: 'text/plain' });
                    const blobUrl = URL.createObjectURL(blob);

                    // Create a link and trigger the download
                    const a = document.createElement('a');
                    a.href = blobUrl;
                    a.download = `${object.id}_${object.seqid}_${seqStart}_${seqEnd}.gb`;//Use the protein id plus the seqid to name the file
                    document.body.appendChild(a);
                    a.click();
                    document.body.removeChild(a);

                    // Revoke the URL
                    URL.revokeObjectURL(blobUrl);
                } catch (error) {
                    console.error('Error fetching data:', error);
                }
            };
            
        function SliderSizes() {
                return (
                    <Box sx={{ width: 150 }}>
                        <Slider
                            size="small"
                            defaultValue={5000}
                            min={0}
                            max={20000}
                            aria-label="Upstream"
                            valueLabelDisplay="auto"
                        />
                        <Slider
                            size="small"
                            defaultValue={5000}
                            min={0}
                            max={20000}
                            aria-label="Downstream"
                            valueLabelDisplay="auto"
                        />
                    </Box>
                );
            }


            return (
                <Fade in={true} timeout={0}> 
                <div 
                    ref={panelRef} 
                    className="info-panel" 
                    style={{ left: `${position.x}px`, top: `${position.y}px` }}
                    onMouseDown={startDrag}
                >
                    <IconButton onClick={onClose} aria-label="close" style={{ position: 'absolute', top: 0, right: 0 }}>
                        <span className="material-icons">close</span>
                    </IconButton>

                    <h3>Information</h3>
                    <p><strong>Protein ID:</strong> {object.id} <a href={`https://www.ncbi.nlm.nih.gov/protein/${object.id}`} target="_blank">🌐</a></p>
                    <p><strong>Cluster:</strong> <span style={{ fontWeight: 'bold', color: `rgb(${object.fillcolor[0]}, ${object.fillcolor[1]}, ${object.fillcolor[2]})` }}>{object.fam_cluster}</span></p>
                    <p><strong>Product:</strong> {object.product}</p>
                        <p><strong>Gene name:</strong> {object.gene}</p>
                        <p><strong>Start:</strong> {object.rel_start}</p>
                        <p><strong>End:</strong> {object.rel_end}</p>
                        <p><strong>Length (aa):</strong> {object.sequence.length}</p>
                        <p><strong>Padloc System:</strong> {object.padloc_system || 'N/A'}</p>
                        <p><strong>Padloc Gene:</strong> {object.padloc_gene || 'N/A'}</p>
                        <p><strong>Deffinder System:</strong> {object.deffinder_system || 'N/A'}</p>
                        <p><strong>Deffinder Gene:</strong> {object.deffinder_gene || 'N/A'}</p>
                        <p><strong>CCTyper System:</strong> {object.cctyper_system || 'N/A'}</p>
                        <p><strong>CCTyper Gene:</strong> {object.cctyper_gene || 'N/A'}</p>
                        <p><strong>Strand:</strong> {object.strand}</p>
                        <p><strong>Nucleotide ID:</strong> {object.seqid} <a href={`https://www.ncbi.nlm.nih.gov/nuccore/${object.seqid}`} target="_blank">🌐</a></p>
                        <p><strong>Assembly ID:</strong> {object.assembly_accession} <a href={`https://www.ncbi.nlm.nih.gov/datasets/genome/${object.assembly_accession}`} target="_blank">🌐</a></p>
                        <p><strong>Species:</strong> {object.species}</p>

                        <pre></pre>


                        <Divider>
                            <Chip label="Actions" size="small" />
                        </Divider>

                        <pre></pre>


                        <div>
                            <div style={{ display: 'flex', justifyContent: 'center' }}>
                                {/* BEGIN: ed8c6549bwf9 */}
                                <Button onClick={() => setShowProteinSequence(!showProteinSequence)}>
                                    {showProteinSequence ? 'Hide' : 'Show'} Protein Sequence
                                </Button>
                                {/* END: ed8c6549bwf9 */}
                            </div>
                            {showProteinSequence && <p>{object.sequence}</p>}
                        </div>

                        <div>
                            <div style={{ display: 'flex', justifyContent: 'center' }}>
                                {/* BEGIN: ed8c6549bwf9 */}
                                <Button  onClick={handleShowStructure}>
                                    {showStructure ? 'Hide' : 'Show'} Protein Structure
                                </Button>
                                {/* END: ed8c6549bwf9 */}
                            </div>
                            {showStructure && <div id="viewport" style={{ width: '300px', height: '300px' }}></div>}
                        </div>

                    {/* Additional Content for Protein Structure and Download Button */}
                    {showStructure && (
                            <div>
                                <div style={{ display: 'flex', justifyContent: 'center' }}>
                                <IconButton onClick={handlePdbDataDownload} aria-label="download">
                                    <span className="material-icons">download</span>
                                </IconButton>
                            </div>

                            </div>
                        )}

                        <div style={{ textAlign: 'center' }}>
                                {showWarning && <p style={{ color: 'red' }}>Proteins above 400aa are not allowed.</p>}
                            </div>

                        {/* New button for showing sliders */}
                        <div style={{ display: 'flex', justifyContent: 'center' }}>
                                {/* BEGIN: ed8c6549bwf9 */}
                        <Button onClick={() => setShowContextSliders(!showContextSliders)}>
                            {showContextSliders ? 'Hide' : 'Get'} Genomic Context
                        </Button>
                                {/* END: ed8c6549bwf9 */}
                            </div>

                        {/* Conditional rendering of sliders */}
                        {showContextSliders && (
                            <Box sx={{ width: 300 }}>
                                <Typography gutterBottom>Upstream</Typography>
                                <Grid container spacing={2} alignItems="center">
                                    <Grid item xs>
                                        <Slider
                                            value={typeof upstreamValue === 'number' ? upstreamValue : 0}
                                            onChange={(e, newValue) => handleSliderChange('upstream', newValue)}
                                            min={0}
                                            max={20000}
                                        />
                                    </Grid>
                                    <Grid item>
                                        <Input
                                            value={upstreamValue}
                                            size="small"
                                            onChange={(e) => handleInputChange(e, 'upstream')}
                                            onBlur={() => handleBlur('upstream')}
                                            inputProps={{
                                                step: 100,
                                                min: 0,
                                                max: 20000,
                                                type: 'number',
                                            }}
                                        />
                                    </Grid>
                                </Grid>
                                <Typography gutterBottom>Downstream</Typography>
                                <Grid container spacing={2} alignItems="center">
                                    <Grid item xs>
                                        <Slider
                                            value={typeof downstreamValue === 'number' ? downstreamValue : 0}
                                            onChange={(e, newValue) => handleSliderChange('downstream', newValue)}
                                            min={0}
                                            max={20000}
                                        />
                                    </Grid>
                                    <Grid item>
                                        <Input
                                            value={downstreamValue}
                                            size="small"
                                            onChange={(e) => handleInputChange(e, 'downstream')}
                                            onBlur={() => handleBlur('downstream')}
                                            inputProps={{
                                                step: 100,
                                                min: 0,
                                                max: 20000,
                                                type: 'number',
                                            }}
                                        />
                                    </Grid>
                                </Grid>
                            </Box>
                        )}

                        {/* Download Button */}

                        <div style={{ display: 'flex', justifyContent: 'center' }}>

                        {showContextSliders && (
                            <IconButton onClick={handleDownload} aria-label="downlaod" >
                            <span className="material-icons">download</span>
                            </IconButton>
                        )}
                        </div>

                    {/* ... Other content ... */}
                </div>
            </Fade>

            );
        }


            ReactDOM.render(<App />, document.getElementById('root'));
        </script>

        </body>
        </html>

        """
        f = open(self.output+"/"+self.output+"_extra.html","w+")
        f.write(first_part+second_part)
        f.close()

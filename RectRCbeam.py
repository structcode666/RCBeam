import streamlit as st
import pandas as pd
import RectRCbeamCalculator as rcbeam
from RectRCbeamCalculator import conc_beam_shear


##Pile Calcualtor Layout##
st.set_page_config(page_title="RC Transfer Rates Calculator", layout="wide")

st.header("RC Beam Rates")
st.markdown("This application helps the user get a feel of the parameters that drive the reinforcement rates of a RC concrete transfer beam")
st.markdown("Our base case will be an banded office floor in a 8.4m x 8.4m grid with a simply supported reinforced concrete transfer band beam. The band beam is transferring a column at mid span supporting 4 office levels. Our base case will be a **700Dx2400W** reinfoced concrete transfer beam that spans **8.0m**")
st.image("RCBeam.png")

with st.sidebar:
    st.header("**Concrete Beam Parameters**")
    beam_depth = st.slider(label="Beam Depth (mm)", min_value = 250, max_value=3000, value=500, step=50, format="%i")
    beam_width = st.slider(label="Beam Width (mm)", min_value = 600, max_value=3600, value=2400, step=600, format="%i")
    beam_span = st.slider(label="Beam Span (m)", min_value = 3, max_value=15, value=8, step=1, format="%f")
    conc_grade = st.selectbox(label="Concrete Grade (MPa)", options=[25, 32, 40, 50, 65, 80], index= 2)

    st.header("Concrete Beam Reinforcement Parameters")

    top_reo = st.selectbox(label="Provide Additional Top reinforcement", options=["Yes", "No"], index= 1)
    ligs = st.selectbox(label="Adjust ligature size and spacing", options=["Yes", "No"], index= 1)

    
    ## Define Conc Beam#
    if top_reo == "Yes" and ligs == "Yes":
        As_top = st.slider(label="Compression Reinforcement (mm2)", min_value = 0, max_value=10000, value=100, step=500, format="%i")
        ## Provide Ligatures##
        lig_size = st.selectbox(label = "Lig Size", options = [10,12,16,20], index = 1)
        lig_sets = st.slider(label="Lig Sets", min_value = 1, max_value=5, value=4, step=1, format="%i")
        lig_spacing = st.slider(label="Lig Spacing (mm)", min_value = 50, max_value=300, value=150, step=30, format="%i")
        conc_beam = conc_beam_shear(b = beam_width, d = beam_depth, conc_grade=conc_grade,lig_size = lig_size, lig_spacing = lig_spacing, lig_legs = 2*lig_sets, fsy_f=500., As_top=As_top )
    elif top_reo == "Yes" and ligs == "No":
        As_top = st.slider(label="Compression Reinforcement  (mm2)", min_value = 0, max_value=10000, value=100, step=500, format="%i")
        conc_beam = conc_beam_shear(b = beam_width, d = beam_depth, conc_grade=conc_grade,lig_size = 12, lig_spacing = 150, lig_legs = 8, fsy_f=500., As_top=As_top)
    elif top_reo == "No" and ligs == "No":
        conc_beam = conc_beam_shear(b = beam_width, d = beam_depth, conc_grade=conc_grade,lig_size = 12, lig_spacing = 150, lig_legs = 8, fsy_f=500.)
    elif top_reo == "No" and ligs == "Yes":
        conc_beam = conc_beam_shear(b = beam_width, d = beam_depth, conc_grade=conc_grade,lig_size = 12, lig_spacing = 150, lig_legs = 8, fsy_f=500.)
        ## Provide Ligatures##
        lig_size = st.selectbox(label = "Lig Size", options = [10,12,16,20], index = 1)
        lig_sets = st.slider(label="Lig Sets", min_value = 1, max_value=5, value=2, step=1, format="%i")
        lig_spacing = st.slider(label="Lig Spacing (mm)", min_value = 50, max_value=300, value=150, step=30, format="%i")
        
        conc_beam = conc_beam_shear(b = beam_width, d = beam_depth, conc_grade=conc_grade,lig_size = lig_size, lig_spacing = lig_spacing, lig_legs = 2*lig_sets, fsy_f=500.)

    st.image("innovislogo.png")

## Define Loads and analyse beam##

load_combo_list = rcbeam.load_combos(1800,850, 55, 25, conc_beam)
conc_beam.des_bending_moment, conc_beam.des_shear, conc_beam.ser_design_moment = rcbeam.beam_analysis(beam_span, load_combo_list)
beam_deflection_dict = rcbeam.beam_deflections(beam_span, conc_beam, load_combo_list)
col1, col2, col3 = st.columns(3)

with col1:
    st.header("Reinforcement Rate")

    if (conc_beam.check_shear() == 'Increase size of beam, or increase ligatures/decrease spacing') or (conc_beam.ult_tensile_reo() == 'Increase beam depth'):
        "Beam depth/reinforcement not adequate"
    else:
        st.metric(label = "", value=f"{rcbeam.calculate_reo_rate(conc_beam)} kg/m3", delta=rcbeam.reo_rate_comparison(conc_beam), delta_color="inverse")

with col2:
    st.header("Deflection") 
    if (conc_beam.check_shear() == 'Increase size of beam, or increase ligatures/decrease spacing') or (conc_beam.ult_tensile_reo() == 'Increase beam depth'):
        "Beam depth/reinforcement not adequate"
    else:
        total_deflection = beam_deflection_dict["Total Deflection"]
        st.metric(label = "", value=f"{total_deflection} mm")
        #  conc_grade
 
 
st.markdown("*This is an educational tool to help design managers get an understanding of the parameters that drive reinforcement rates of a reinforced concrete transfer beam. Reinforcement rates and deflection values have been provided as a guide only and must not be used for final design. Please contact Innovis if project specific designs are required.*")




























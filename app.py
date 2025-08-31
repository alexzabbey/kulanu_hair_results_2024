import streamlit as st
import pandas as pd
import geopandas as gpd
import pydeck as pdk
import altair as alt
import matplotlib.colors as mcolors
from io import BytesIO
import numpy as np


# polygon_layer = pdk.Layer(
#     "PolygonLayer",
#     LAND_COVER,
#     stroked=False,
#     # processes the data as a flat longitude-latitude pair
#     get_polygon="-",
#     get_fill_color=[0, 0, 0, 20],
# )

st.markdown(
    """
<style>
body, html {
    direction: rtl;
    text-align: right;
}
p, div, input, label, h1, h2, h3, h4, h5, h6 {
    direction: rtl;
    text-align: right;
}
h1 > div, h2 > div,h3 > div,  h4, h5, h6, p{
    left: 0 !important;
}
</style>
""",
    unsafe_allow_html=True,
)


def absolute_map(display):
    if nei := ("שכונות" in display):
        if(display == "שכונות קול העיר"):
            gdf = pd.DataFrame(gpd.read_file("kol_results_2024_by_neighborhood.geojson"))
        else:
            year = display.split(" ")[1]
            gdf = pd.DataFrame(gpd.read_file(f"results_{year}_by_neighborhood.geojson"))
        
        gdf["fill"] = gdf["fill"].apply(lambda x: [x * 255 for x in mcolors.to_rgb(x)])
        gdf["geometry"] = gdf["geometry"].apply(lambda x: list(x.exterior.coords))
        layer = pdk.Layer(
            "PolygonLayer",
            gdf,
            id="geojson",
            opacity=0.7,
            stroked=True,
            get_polygon="geometry",
            filled=True,
            get_fill_color="fill",
            get_line_color=[200, 200, 200],
            line_width_min_pixels=2,
            auto_highlight=True,
            pickable=True,
        )
    elif display == "אשכולות 2024":
        eshkolot = gpd.read_file("results_2024_by_eshkol.geojson")
        eshkolot["longitude"] = eshkolot["geometry"].apply(lambda x: x.coords[0][0])
        eshkolot["latitude"] = eshkolot["geometry"].apply(lambda x: x.coords[0][1])
        layer = pdk.Layer(
            "ScatterplotLayer",
            eshkolot,
            get_position=["longitude", "latitude"],
            get_fill_color=[255, 0, 0],
            get_radius=100,
            pickable=True,
            auto_highlight=True,
        )
    elif display == "אשכולות קול העיר":
        eshkolot = gpd.read_file("kol_results_2024_by_eshkol.geojson")
        eshkolot["longitude"] = eshkolot["geometry"].apply(lambda x: x.coords[0][0])
        eshkolot["latitude"] = eshkolot["geometry"].apply(lambda x: x.coords[0][1])
        layer = pdk.Layer(
            "ScatterplotLayer",
            eshkolot,
            get_position=["longitude", "latitude"],
            get_fill_color=[255, 0, 0],
            get_radius=100,
            pickable=True,
            auto_highlight=True,
        )

    datum = "neighborhood" if nei else "name"
    field = "שכונה" if nei else "אשכול"

    tooltip = {"html": f"<b>{field}:</b> {{{datum}}} <br /> <b>קולות:</b> {{ain}}"}
    view_state = pdk.ViewState(
        **{
            "latitude": 32.084759,
            "longitude": 34.781635,
            "zoom": 11.3,
            "maxZoom": 16,
            "bearing": 0,
        }
    )
    return pdk.Deck(
        layer,
        initial_view_state=view_state,
        map_style=pdk.map_styles.LIGHT,
        tooltip=tooltip,
    )

def party_map(filename, display_name, value_column=None):
    # Load the GeoJSON
    gdf = gpd.read_file(filename)
    gdf["longitude"] = gdf["geometry"].apply(lambda x: x.coords[0][0])
    gdf["latitude"] = gdf["geometry"].apply(lambda x: x.coords[0][1])

    # Auto-detect SUM_* column if not given
    if value_column is None:
        sum_cols = [c for c in gdf.columns if c.startswith("SUM_")]
        if sum_cols:
            value_column = sum_cols[0]
        else:
            value_column = "COUNT_kalpi_number"  # fallback

    # Normalize vote counts
    min_votes = gdf[value_column].min()
    max_votes = gdf[value_column].max()

    # RGB for gradient
    rgb_light_green = np.array([144, 238, 144])
    rgb_dark_blue = np.array([0, 0, 139])

    def vote_to_gradient(v):
        # Linear interpolation
        ratio = (v - min_votes) / (max_votes - min_votes) if max_votes > min_votes else 0
        color = rgb_light_green + ratio * (rgb_dark_blue - rgb_light_green)
        return [int(c) for c in color]

    gdf["color"] = gdf[value_column].apply(vote_to_gradient)

    # ScatterplotLayer
    layer = pdk.Layer(
        "ScatterplotLayer",
        gdf,
        get_position=["longitude", "latitude"],
        get_fill_color="color",
        get_radius=120,
        pickable=True,
        auto_highlight=True,
    )

    view_state = pdk.ViewState(
        latitude=32.084759,
        longitude=34.781635,
        zoom=11.3,
        maxZoom=16,
        bearing=0,
    )

    tooltip = {"html": f"<b>קלפי:</b> {{name}} <br/> <b>{display_name}:</b> {{{value_column}}}"}

    return pdk.Deck(
        layer,
        initial_view_state=view_state,
        map_style=pdk.map_styles.LIGHT,
        tooltip=tooltip,
    )



parties = [
    #{"file": "results_all.json", "name": "כל הרשימות", "vote_col": "SUM_all"},
    {"file": "results_ain.json", "name": "כולנו העיר", "vote_col": "SUM_ain"},
    {"file": "results_t.json", "name": "תקווה חדשה", "vote_col": "SUM_t"},
    {"file": "results_shas.json", "name": "ש״ס", "vote_col": "SUM_shas"},
    {"file": "results_likud.json", "name": "הליכוד", "vote_col": "SUM_likud"},
    {"file": "results_ph.json", "name": "יש עתיד", "vote_col": "SUM_ph"},
    {"file": "results_meretz.json", "name": "מרצ", "vote_col": "SUM_meretz"},
    {"file": "results_ta.json", "name": "תל אביב 1", "vote_col": "SUM_ta"},
    {"file": "results_df.json", "name": "תושבי העיר", "vote_col": "SUM_df"},
    {"file": "results_dn.json", "name": "הרשימה השווה", "vote_col": "SUM_dn"},
    {"file": "results_gb.json", "name": "מאמינים בתל אביב", "vote_col": "SUM_gb"},
    {"file": "results_hi.json", "name": "חי - חילונים ירוקים", "vote_col": "SUM_hi"},
    {"file": "results_kl.json", "name": "קול העיר", "vote_col": "SUM_kl"},
    {"file": "results_r.json", "name": "חרות", "vote_col": "SUM_r"},
    {"file": "results_tz.json", "name": "סיעת הצעירים", "vote_col": "SUM_tz"},
    {"file": "results_zk.json", "name": "סיעת הגימלאים והאזרחים הוותיקים", "vote_col": "SUM_zk"}
]


def neighborhood_party_map(joined_file, display_name, vote_col):
    gdf = gpd.read_file(joined_file)
    
    if gdf.crs.is_geographic:
        gdf = gdf.to_crs(epsg=3857)

    if gdf.crs is None:
        gdf = gdf.set_crs(epsg=3857)  # your file is in Web Mercator
    gdf = gdf.to_crs(epsg=4326)
    
    min_votes = gdf[vote_col].min()
    max_votes = gdf[vote_col].max()
    
    def vote_to_color(v):
        if v is None or pd.isna(v):
            return [0, 0, 0, 0]
        ratio = (v - min_votes) / (max_votes - min_votes) if max_votes > min_votes else 0
        r = int(144 * (1 - ratio) + 0 * ratio)
        g = int(238 * (1 - ratio) + 0 * ratio)
        b = int(144 * (1 - ratio) + 139 * ratio)
        return [r, g, b]

    gdf["fill"] = gdf[vote_col].apply(vote_to_color)
    
    def poly_coords(geom):
        if geom.type == "Polygon":
            return [list(geom.exterior.coords)]
        elif geom.type == "MultiPolygon":
            return [list(p.exterior.coords) for p in geom.geoms]
        return []
    
    gdf["geometry_coords"] = gdf["geometry"].apply(poly_coords)
    
    layer = pdk.Layer(
        "PolygonLayer",
        gdf,
        get_polygon="geometry_coords",
        get_fill_color="fill",
        get_line_color=[200, 200, 200],
        stroked=True,
        filled=True,
        pickable=True,
        auto_highlight=True,
        opacity=0.7,
        line_width_min_pixels=2,
    )

    tooltip = {"html": f"<b>שכונה:</b> {{neighborhood}} <br/> <b>{display_name}:</b> {{{vote_col}}}"}
    view_state = pdk.ViewState(latitude=32.084759, longitude=34.781635, zoom=11.3)

    return pdk.Deck(layers=[layer], initial_view_state=view_state, map_style=pdk.map_styles.LIGHT, tooltip=tooltip)




def diff():
    joined_final = gpd.read_file("results_2024_by_neighborhood.geojson")
    return (
        alt.Chart(joined_final.drop(columns=["geometry"]))
        .mark_bar()
        .encode(
            x=alt.X(
                "diff:Q",
                title="השינוי בין 2018 ל-2024 באחוז הקולות מתוך כלל הקולות שהצביעו לנו בשכונה זו",
                axis=alt.Axis(format=".0%", tickCount=7),
            ),
            y=alt.Y("neighborhood:N", title="Neighborhood", sort="-x"),
            color=alt.Color("diff:Q", scale=alt.Scale(scheme="redblue")),
            tooltip=["neighborhood", "diff"],
        )
        .interactive()
    )  # .properties(width=800, height=800)


def full():
    gdf = gpd.read_file("results_2024_and_2018.geojson")
    gdf["כמות קולות לכולנו העיר 2018"] = gdf["כמות קולות לכולנו העיר 2018"].astype(int)
    gdf["כמות קולות 2018"] = gdf["כמות קולות 2018"].astype(int)
    gdf["בעלי זכות בחירה 2018"] = gdf["בעלי זכות בחירה 2018"].astype(int)
    return gdf


st.write("## ניתוח תוצאות הבחירות המקומיות 2024 - כולנו העיר")
st.write(
    "היי, עשיתי כמה ניתוחים של נתוני הבחירות. התוצאות של 2024 מדויקות, אבל ייתכנו טעויות בנתונים של 2018, שעברתי עליהם בצורה קצת פחות מדוקדקת. אם יש עוד ניתוחים שמעניין אתכם לראות מוזמנים לכתוב לי :). -- אלכס"
)
st.write("### כמות קולות אבסולוטית לפי שכונות")
st.write(
    "מפה של כמות הקולות האבסולוטית בכל שכונה. שכונות שלא מופיעות הן כאלה שלא היו בהן אשכולות קלפיות"
)
year = st.radio(
    "בחרו", ["שכונות 2024", "אשכולות 2024", "שכונות 2018", "אשכולות קול העיר", "שכונות קול העיר"], index=0, horizontal=True
)
st.pydeck_chart(absolute_map(year))
st.write("### מפות לפי רשימות - 2024")

# Dropdown for neighborhood map
party_names = [p["name"] for p in parties]
selected_name = st.selectbox("בחרו רשימה למפה לפי שכונות:", [p["name"] for p in parties])
selected_party = next(p for p in parties if p["name"] == selected_name)

st.pydeck_chart(neighborhood_party_map(
    "all_parties_by_neighborhood.geojson",
    selected_party["name"],
    selected_party["vote_col"]   
))

st.markdown(
    """
    <div style="display: flex; align-items: center; gap: 10px; margin-top: 10px;">
        <span>High votes</span>
        <div style="
            width: 200px;
            height: 20px;
            background: linear-gradient(to right, rgb(144,238,144), rgb(0,0,139));
            border: 1px solid #000;
        "></div>
        <span>Low votes</span>
    </div>
    """,
    unsafe_allow_html=True
)

st.write("")
st.write("")

# Build dropdown options
selected_name = st.selectbox("בחרו רשימה למפה לפי אשכולות:", [p["name"] for p in parties])
selected_party = next(p for p in parties if p["name"] == selected_name)

# Find the selected party info
selected_party = next(p for p in parties if p["name"] == selected_name)

# Display the map
st.pydeck_chart(
    party_map(
        selected_party["file"],
        selected_party["name"],
        value_column=selected_party["vote_col"]
    )
)


st.write("---")
st.write("###  שינוי באחוז הקולות לפי שכונות")

st.write(
    "הגרף מראה את השינוי באחוז הקולות מתוך כלל הקולות שהצביעו לנו בין הבחירות האחרונות לבחירות הקודמות. כלומר, בכחול שגונות שהתחזקו ובאדום שכונות שנחלשו."
)
st.altair_chart(diff(), use_container_width=True)
st.write("---")
with st.expander("התוצאות המלאות", expanded=False):
    gdf = pd.DataFrame(full())
    gdf["geometry"] = gdf["geometry"].astype(str)
    # st.download_button("להורדת הנתונים בפורמט geojson", gdf.to_file(), "results.geojson"
    excel_file = BytesIO()
    gdf.to_excel(excel_file, index=False)
    excel_file.seek(0)
    st.download_button(
        "להורדת הנתונים באקסל",
        excel_file,
        "results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
    st.table(gdf)


# לסמן על המפה איפה הצבנו אנשים
# לבדוק את קולות לפי שכונות 2018

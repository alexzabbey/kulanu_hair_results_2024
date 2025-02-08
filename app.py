import streamlit as st
import pandas as pd
import geopandas as gpd
import pydeck as pdk
import altair as alt
import matplotlib.colors as mcolors
from io import BytesIO


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
st.write("---")
st.write("### שינוי באחוז הקולות לפי שכונות")

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

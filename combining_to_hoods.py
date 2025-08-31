import geopandas as gpd
import os

# Folder containing your Eshkol-level JSON files
eshkol_folder = r"C:\Users\ShiraTomer_Laptop\Documents\standingTogether\kulanu_hair_results_2024"

# Neighborhood polygons
neighborhood_file = "neighborhoods.geojson"

# Output file
output_file = "all_parties_by_neighborhood.geojson"

# List of party files and vote columns
parties = [
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

# Load neighborhoods
neighborhoods = gpd.read_file(neighborhood_file)
if neighborhoods.crs.is_geographic:
    neighborhoods = neighborhoods.to_crs(epsg=3857)

# Start with just the neighborhoods GeoDataFrame
joined_gdf = neighborhoods.copy()

for party in parties:
    print(f"Processing {party['file']}...")
    eshkol = gpd.read_file(os.path.join(eshkol_folder, party["file"]))
    
    # Use centroids if not points
    if eshkol.geometry.geom_type.iloc[0] != "Point":
        eshkol["centroid"] = eshkol.geometry.centroid
        eshkol_points = eshkol.set_geometry("centroid")
    else:
        eshkol_points = eshkol

    # Ensure same CRS
    if eshkol_points.crs != joined_gdf.crs:
        eshkol_points = eshkol_points.to_crs(joined_gdf.crs)

    # Spatial join
    joined_temp = gpd.sjoin(eshkol_points, neighborhoods, how="left", predicate="within")

    # Aggregate votes by neighborhood
    neighborhood_votes = joined_temp.groupby("neighborhood")[party["vote_col"]].sum().reset_index()

    # Merge into master GeoDataFrame
    joined_gdf = joined_gdf.merge(neighborhood_votes, on="neighborhood", how="left")

# Save the master file
joined_gdf.to_file(output_file, driver="GeoJSON")
print(f"Saved {output_file}")

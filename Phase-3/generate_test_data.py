# # #!/usr/bin/env python3

# # import json

# # import math

# # import random

# # import argparse



# # # ----------------- Geo helpers -----------------



# # def haversine_m(lat1, lon1, lat2, lon2):

# #     """

# #     Haversine distance in meters, lat/lon in degrees.

# #     """

# #     R = 6371000.0

# #     phi1 = math.radians(lat1)

# #     phi2 = math.radians(lat2)

# #     dphi = math.radians(lat2 - lat1)

# #     dlambda = math.radians(lon2 - lon1)



# #     a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2

# #     c = 2 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))

# #     return R * c



# # # ----------------- Graph generation -----------------



# # def generate_random_nodes(num_nodes, base_lat, base_lon, radius_km):

# #     """

# #     Generate nodes scattered uniformly within a circle of radius_km

# #     around (base_lat, base_lon). This looks like a "city" area.

# #     """

# #     nodes = []

# #     for i in range(num_nodes):

# #         # Random point in circle (radius in km)

# #         r = radius_km * math.sqrt(random.random())

# #         theta = random.uniform(0, 2 * math.pi)



# #         # Convert radial offset in km to degrees:

# #         # approx 1 deg lat ~ 111 km; 1 deg lon ~ 111 km * cos(lat)

# #         dlat = (r * math.cos(theta)) / 111.0

# #         dlon = (r * math.sin(theta)) / (111.0 * math.cos(math.radians(base_lat)))



# #         lat = base_lat + dlat

# #         lon = base_lon + dlon



# #         nodes.append({"id": i, "lat": lat, "lon": lon})

# #     return nodes





# # def build_random_geometric_graph(num_nodes, base_lat, base_lon, radius_km,

# #                                  k_nearest=6, min_speed_kmh=20.0, max_speed_kmh=60.0):

# #     """

# #     Build a "realistic-ish" road graph:

# #       - nodes randomly scattered in a city-radius

# #       - each node connects to k_nearest neighbors (undirected)

# #       - edge length from haversine, average_time from random speed

# #     """

# #     nodes = generate_random_nodes(num_nodes, base_lat, base_lon, radius_km)



# #     # Precompute positions

# #     positions = {n["id"]: (n["lat"], n["lon"]) for n in nodes}



# #     # For each node, find k nearest neighbors

# #     adjacency = {i: set() for i in range(num_nodes)}

# #     for i in range(num_nodes):

# #         lat_i, lon_i = positions[i]

# #         dists = []

# #         for j in range(num_nodes):

# #             if i == j:

# #                 continue

# #             lat_j, lon_j = positions[j]

# #             d = haversine_m(lat_i, lon_i, lat_j, lon_j)

# #             dists.append((d, j))

# #         dists.sort(key=lambda x: x[0])

# #         for _, j in dists[:k_nearest]:

# #             adjacency[i].add(j)

# #             adjacency[j].add(i)  # undirected



# #     # Now create edge list

# #     edges = []

# #     edge_id = 0

# #     seen = set()

# #     for u in range(num_nodes):

# #         for v in adjacency[u]:

# #             if (v, u) in seen:

# #                 continue

# #             lat_u, lon_u = positions[u]

# #             lat_v, lon_v = positions[v]

# #             length = haversine_m(lat_u, lon_u, lat_v, lon_v)  # meters



# #             # speed in km/h -> m/s

# #             speed_kmh = random.uniform(min_speed_kmh, max_speed_kmh)

# #             speed_mps = speed_kmh * 1000.0 / 3600.0

# #             avg_time = max(length / speed_mps, 1.0)  # seconds



# #             edges.append({

# #                 "id": edge_id,

# #                 "u": u,

# #                 "v": v,

# #                 "length": length,

# #                 "average_time": avg_time,

# #                 "oneway": False,

# #                 "type": "road"

# #             })

# #             edge_id += 1

# #             seen.add((u, v))



# #     meta = {

# #         "nodes": num_nodes,

# #         "edges": len(edges)

# #     }



# #     return {

# #         "meta": meta,

# #         "nodes": nodes,

# #         "edges": edges

# #     }



# # # ----------------- Queries generation -----------------



# # def generate_queries(num_nodes: int,

# #                      num_events: int,

# #                      min_orders: int,

# #                      max_orders: int,

# #                      max_drivers: int):

# #     """

# #     Generate queries.json with extra fields:

# #       order_price, is_cancelled, is_premium

# #     """

# #     events = []

# #     for eid in range(num_events):

# #         num_orders = random.randint(min_orders, max_orders)

# #         num_drivers = random.randint(1, max_drivers)

# #         depot_node = random.randrange(num_nodes)



# #         orders = []

# #         for oid in range(num_orders):

# #             pickup = random.randrange(num_nodes)

# #             dropoff = random.randrange(num_nodes)

# #             while dropoff == pickup:

# #                 dropoff = random.randrange(num_nodes)



# #             price = round(random.uniform(100.0, 1000.0), 2)

# #             is_premium = random.random() < 0.25    # 25% premium

# #             is_cancelled = random.random() < 0.05  # 5% cancelled



# #             orders.append({

# #                 "order_id": oid,

# #                 "pickup": pickup,

# #                 "dropoff": dropoff,

# #                 "order_price": price,

# #                 "is_cancelled": is_cancelled,

# #                 "is_premium": is_premium

# #             })



# #         events.append({

# #             "event_id": eid,

# #             "orders": orders,

# #             "fleet": {

# #                 "num_delievery_guys": num_drivers,  # matches your C++ key

# #                 "depot_node": depot_node

# #             }

# #         })



# #     return {

# #         "meta": {

# #             "num_events": num_events

# #         },

# #         "events": events

# #     }



# # # ----------------- Main -----------------



# # def main():

# #     parser = argparse.ArgumentParser(description="Generate realistic-ish geo graph + queries.")

# #     parser.add_argument("--nodes", type=int, default=100,

# #                         help="Number of graph nodes.")

# #     parser.add_argument("--events", type=int, default=10,

# #                         help="Number of events in queries.json.")

# #     parser.add_argument("--min-orders", type=int, default=5,

# #                         help="Minimum number of orders per event.")

# #     parser.add_argument("--max-orders", type=int, default=30,

# #                         help="Maximum number of orders per event.")

# #     parser.add_argument("--max-drivers", type=int, default=10,

# #                         help="Maximum number of drivers per event.")

# #     parser.add_argument("--seed", type=int, default=42,

# #                         help="Random seed for reproducibility.")



# #     # City-like region parameters

# #     parser.add_argument("--base-lat", type=float, default=19.07,

# #                         help="Center latitude (e.g., 19.07 ~ Mumbai).")

# #     parser.add_argument("--base-lon", type=float, default=72.88,

# #                         help="Center longitude (e.g., 72.88 ~ Mumbai).")

# #     parser.add_argument("--radius-km", type=float, default=10.0,

# #                         help="Radius of area in km (city size).")



# #     parser.add_argument("--graph-out", type=str, default="graph.json",

# #                         help="Output file for graph JSON.")

# #     parser.add_argument("--queries-out", type=str, default="queries.json",

# #                         help="Output file for queries JSON.")



# #     args = parser.parse_args()

# #     random.seed(args.seed)



# #     print(f"[*] Generating graph with {args.nodes} nodes "

# #           f"around ({args.base_lat}, {args.base_lon}) radius {args.radius_km} km...")

# #     graph = build_random_geometric_graph(

# #         num_nodes=args.nodes,

# #         base_lat=args.base_lat,

# #         base_lon=args.base_lon,

# #         radius_km=args.radius_km,

# #         k_nearest=6

# #     )



# #     print(f"[*] Generating {args.events} events with orders and extra fields...")

# #     queries = generate_queries(

# #         num_nodes=args.nodes,

# #         num_events=args.events,

# #         min_orders=args.min_orders,

# #         max_orders=args.max_orders,

# #         max_drivers=args.max_drivers

# #     )



# #     with open(args.graph_out, "w") as f:

# #         json.dump(graph, f, indent=2)



# #     with open(args.queries_out, "w") as f:

# #         json.dump(queries, f, indent=2)



# #     print(f"[+] Wrote {args.graph_out} and {args.queries_out}")





# # if __name__ == "__main__":

# #     main()










































# import json
# import math
# import random
# import copy
# from typing import List, Dict, Any, Tuple

# # --- Configuration Constants (From C++ scheduler.cpp) ---
# SEED = 42
# MAX_PRICE = 1000.0
# PREMIUM_BENEFIT = 5
# PRICE_NORMALIZER = 10
# K_MEANS_ITERATIONS = 10
# EPSILON_DISTANCE = 1e-5
# DEPOT_NODE = 0 # Default depot node from C++ init


# # --- 1. GEO HELPERS (From original Python generator) ---
# def haversine_m(lat1, lon1, lat2, lon2):
#     """Haversine distance in meters, lat/lon in degrees."""
#     R = 6371000.0
#     phi1 = math.radians(lat1)
#     phi2 = math.radians(lat2)
#     dphi = math.radians(lat2 - lat1)
#     dlambda = math.radians(lon2 - lon1)

#     a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
#     c = 2 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
#     return R * c

# # --- 2. C++ STRUCTS / CLASSES IN PYTHON ---

# class Order:
#     def __init__(self, data: Dict[str, Any]):
#         self.id = data["order_id"]
#         self.pickupNode = data["pickup"]
#         self.dropNode = data["dropoff"]
#         self.order_price = data["order_price"]
#         self.is_premium = data["is_premium"]
#         self.is_cancelled = data["is_cancelled"]
#         self.priority_score = 0
#         self.delivered = False
#         self.completionTime = 0.0

# class Graph:
#     """A minimal Graph class to provide node locations and euclidean distance."""
#     def __init__(self, nodes_data: List[Dict[str, Any]]):
#         self.nodes = {n['id']: n for n in nodes_data}
#         # Pre-calculating a Max Speed proxy for travel time estimation
#         self.MAX_SPEED_MPS = 60.0 * 1000.0 / 3600.0 # 60km/h in m/s

#     def euclidean_estimate(self, u: int, v: int) -> float:
#         """Euclidean distance in meters, used as a fast heuristic."""
#         if u not in self.nodes or v not in self.nodes:
#             return 0.0
#         n_u = self.nodes[u]
#         n_v = self.nodes[v]
#         return haversine_m(n_u['lat'], n_u['lon'], n_v['lat'], n_v['lon'])

#     def travel_time_proxy(self, u: int, v: int) -> float:
#         """Simple distance / max_speed proxy for A* estimation."""
#         distance = self.euclidean_estimate(u, v)
#         # Add a minimum time buffer (similar to EPSILON in C++ or a simple delay)
#         return max(distance / self.MAX_SPEED_MPS, 5.0)

# class Driver:
#     def __init__(self, driver_id: int, depot_node: int):
#         self.driver_id = driver_id
#         self.currentNode = depot_node
#         self.currentTime = 0.0
#         self.route = [depot_node]
#         self.actual_path = [depot_node] # Simplified to match route for this proxy
#         self.order_ids = []
#         self.carryingBag = []
#         self.total_distance_m = 0.0

# # --- 3. C++ SCHEDULER LOGIC IMPLEMENTATION (in Python) ---

# def calculate_priority(orders: List[Order]):
#     """Implements the C++ priority calculation logic."""
#     for ord in orders:
#         if ord.is_cancelled:
#             ord.priority_score = 0
#             continue
        
#         # a = floor((price / max_price) * 10)
#         a = math.floor((ord.order_price / MAX_PRICE) * PRICE_NORMALIZER)
#         # b = 5 (premium) or 1 (normal)
#         b = PREMIUM_BENEFIT if ord.is_premium else 1
#         ord.priority_score = a * b

# def kmeans_clustering(K: int, orders: List[Order], graph: Graph) -> List[List[Order]]:
#     """Implements the C++ K-Means clustering logic to assign orders."""
#     valid_orders = [o for o in orders if not o.is_cancelled]
#     if not valid_orders:
#         return [[] for _ in range(K)]

#     # If fewer orders than drivers, assign one order per driver (simplification)
#     if len(valid_orders) <= K:
#         assignments = [[] for _ in range(K)]
#         for i, order in enumerate(valid_orders):
#             assignments[i].append(order)
#         return assignments

#     def get_midpoint(order: Order) -> Tuple[float, float]:
#         """Calculates the midpoint lat/lon for an order."""
#         p_node = graph.nodes[order.pickupNode]
#         d_node = graph.nodes[order.dropNode]
#         return (p_node['lat'] + d_node['lat']) / 2.0, (p_node['lon'] + d_node['lon']) / 2.0

#     def distance_sq(lat1, lon1, lat2, lon2):
#         """Euclidean distance squared in the lat/lon plane."""
#         return (lat1 - lat2) ** 2 + (lon1 - lon2) ** 2

#     # Initialization: K-means++ style (as implemented in C++)
#     centroids: List[Tuple[float, float]] = []
    
#     # Start with the first order's midpoint
#     centroids.append(get_midpoint(valid_orders[0])) 
    
#     for i in range(1, K):
#         max_dist = -1.0
#         best_candidate_order: Order = valid_orders[0]
        
#         # Find the order whose midpoint is furthest from all currently chosen centroids
#         for order in valid_orders:
#             oLat, oLon = get_midpoint(order)
#             min_dist_to_centroid = float('inf')
            
#             for j in range(len(centroids)):
#                 cLat, cLon = centroids[j]
#                 d = distance_sq(oLat, oLon, cLat, cLon)
#                 min_dist_to_centroid = min(min_dist_to_centroid, d)
            
#             if min_dist_to_centroid > max_dist:
#                 max_dist = min_dist_to_centroid
#                 best_candidate_order = order
        
#         centroids.append(get_midpoint(best_candidate_order))


#     # K-Means Iterations
#     for i in range(K_MEANS_ITERATIONS):
#         assignments: List[List[Order]] = [[] for _ in range(K)]
#         new_centroids: List[Tuple[float, float]] = [(0.0, 0.0)] * K
#         counts = [0] * K

#         # Assignment Step
#         for order in valid_orders:
#             oLat, oLon = get_midpoint(order)
#             bestDist = float('inf')
#             bestK = 0

#             for k in range(K):
#                 cLat, cLon = centroids[k]
#                 d = distance_sq(oLat, oLon, cLat, cLon)
#                 if d < bestDist:
#                     bestDist = d
#                     bestK = k
            
#             assignments[bestK].append(order)
#             new_centroids[bestK] = (new_centroids[bestK][0] + oLat, new_centroids[bestK][1] + oLon)
#             counts[bestK] += 1

#         # Update Step
#         changed = False
#         for k in range(K):
#             if counts[k] > 0:
#                 new_lat = new_centroids[k][0] / counts[k]
#                 new_lon = new_centroids[k][1] / counts[k]
                
#                 old_lat, old_lon = centroids[k]
#                 if abs(new_lat - old_lat) > EPSILON_DISTANCE or abs(new_lon - old_lon) > EPSILON_DISTANCE:
#                     changed = True
                
#                 centroids[k] = (new_lat, new_lon)
        
#         if not changed:
#             break
            
#     return assignments

# def solve_driver_route_dummy(driver: Driver, assigned_orders: List[Order], graph: Graph):
#     """
#     A simplified VRP simulation to generate route and metrics.
    
#     Logic: Sort orders by priority score (descending) and create a route 
#     by greedily adding (Pickup, Dropoff) pairs.
#     """
    
#     # 1. Sort orders by priority (highest first)
#     assigned_orders.sort(key=lambda o: o.priority_score, reverse=True)
    
#     route = [driver.currentNode]
#     current_node = driver.currentNode
#     current_time = 0.0
#     total_distance = 0.0
    
#     for order in assigned_orders:
#         if order.is_cancelled:
#             continue

#         # Add pickup
#         total_distance += graph.euclidean_estimate(current_node, order.pickupNode)
#         current_time += graph.travel_time_proxy(current_node, order.pickupNode)
#         route.append(order.pickupNode)
#         current_node = order.pickupNode
#         driver.carryingBag.append(order.id) # Pickup event

#         # Add dropoff
#         total_distance += graph.euclidean_estimate(current_node, order.dropNode)
#         current_time += graph.travel_time_proxy(current_node, order.dropNode)
#         route.append(order.dropNode)
#         current_node = order.dropNode
        
#         # Mark as delivered and calculate completion time
#         order.delivered = True
#         order.completionTime = current_time
#         driver.order_ids.append(order.id)
#         driver.carryingBag.remove(order.id) # Dropoff event

#     # Return to depot
#     total_distance += graph.euclidean_estimate(current_node, driver.currentNode)
#     current_time += graph.travel_time_proxy(current_node, driver.currentNode)
#     route.append(driver.currentNode)
    
#     driver.route = route
#     driver.actual_path = route # Simplified: actual_path = route for proxy
#     driver.currentTime = current_time
#     driver.total_distance_m = total_distance

# def process_query_in_python(event_data: Dict[str, Any], graph: Graph) -> Dict[str, Any]:
#     """Simulates the C++ Scheduler::process_query function."""
    
#     # --- 1. Order Processing and Priority Calculation ---
#     orders: List[Order] = []
#     for o_data in event_data["orders"]:
#         orders.append(Order(o_data))
    
#     calculate_priority(orders)
    
#     # --- 2. Fleet Parsing and Initialization ---
#     numDrivers = event_data["fleet"]["num_delievery_guys"]
#     depotNode = event_data["fleet"]["depot_node"]
#     drivers: List[Driver] = [Driver(i, depotNode) for i in range(numDrivers)]
    
#     # --- 3. K-Means Clustering ---
#     assigned_orders_groups: List[List[Order]] = kmeans_clustering(numDrivers, orders, graph)
    
#     # --- 4. Solve Driver Routes (Dummy VRP) ---
#     for i in range(numDrivers):
#         if i < len(assigned_orders_groups) and assigned_orders_groups[i]:
#             solve_driver_route_dummy(drivers[i], assigned_orders_groups[i], graph)
#         else:
#             # Driver assigned no orders just stays at depot
#             drivers[i].route = [depotNode]
#             drivers[i].actual_path = [depotNode]

#     # --- 5. Metrics Calculation and Output Formatting ---
    
#     totalDeliveryTime = sum(o.completionTime for o in orders if o.delivered)
#     total_orders_delivered = sum(1 for o in orders if o.delivered)
    
#     # Re-calculate required metrics for the structure specified in previous turn
#     total_revenue = sum(o.order_price for o in orders if o.delivered)
#     total_cancelled_orders = sum(1 for o in orders if o.is_cancelled)
#     average_delivery_time = totalDeliveryTime / total_orders_delivered if total_orders_delivered > 0 else 0.0
    
#     output_drivers = []
#     for d in drivers:
#         if d.order_ids:
#             # Note: total_distance_m is an internal proxy for better dummy data
#             output_drivers.append({
#                 "driver_id": d.driver_id,
#                 "route_nodes": d.route, # This is Driver::route/actual_path from C++
#                 "orders_delivered": d.order_ids,
#                 "total_travel_time": round(d.currentTime, 3), 
#                 "total_distance_m": round(d.total_distance_m, 2)
#             })

#     # The C++ main.cpp structure includes "assignments" and "metrics"
#     # The structure suggested in the previous turn merges this for clarity (here we stick to the output format in C++)

#     assignments = []
#     for d in drivers:
#          assignments.append({
#             "driver_id": d.driver_id,
#             "route": d.route,
#             "order_ids": d.order_ids
#          })

#     # This mirrors the C++ main.cpp and scheduler.cpp structure
#     result = {
#         "event_id": event_data["event_id"],
#         # processing_time is added by main.cpp (simulated later)
#         "assignments": assignments,
#         "metrics": {
#             "total_delivery_time_s": round(totalDeliveryTime, 3), # Matches C++ metric key
#             "total_orders_delivered": total_orders_delivered, # Added for richer output structure
#             "total_revenue": round(total_revenue, 2), # Added for richer output structure
#             "total_cancelled_orders": total_cancelled_orders, # Added for richer output structure
#             "average_delivery_time": round(average_delivery_time, 3) # Added for richer output structure
#         }
#     }
    
#     # Add the drivers structure we had in the previous turn for completeness/richer output
#     # Note: C++ output for "results" is `result` (which has assignments/metrics). 
#     # I will add the detailed driver metrics back into the main structure for completeness.
#     # I will adjust the final output to match the C++ structure but use the richer driver info.

#     final_result_for_output = {
#         "event_id": event_data["event_id"],
#         # processing_time is added by main.cpp (simulated later)
#         "assignments": assignments, # The C++ output structure
#         "drivers": output_drivers, # Added for rich output validation
#         "metrics": {
#             "total_delivery_time_s": round(totalDeliveryTime, 3),
#             "total_orders_delivered": total_orders_delivered,
#             "total_revenue": round(total_revenue, 2),
#             "total_cancelled_orders": total_cancelled_orders,
#             "average_delivery_time": round(average_delivery_time, 3)
#         }
#     }

#     return final_result_for_output

# def generate_expected_output_json():
#     """Generates graph.json and queries.json inputs, then processes them."""
#     random.seed(SEED)
    
#     # --- Generate Input (Graph & Queries) using original Python logic ---
    
#     # Default arguments from the user's main block
#     num_nodes = 100
#     num_events = 10
#     min_orders = 5
#     max_orders = 30
#     max_drivers = 10
#     base_lat = 19.07
#     base_lon = 72.88
#     radius_km = 10.0
#     k_nearest = 6

#     # 1. Graph Generation (Minimal for node location)
#     def generate_random_nodes(num_nodes, base_lat, base_lon, radius_km):
#         nodes = []
#         for i in range(num_nodes):
#             r = radius_km * math.sqrt(random.random())
#             theta = random.uniform(0, 2 * math.pi)
#             dlat = (r * math.cos(theta)) / 111.0
#             dlon = (r * math.sin(theta)) / (111.0 * math.cos(math.radians(base_lat)))
#             lat = base_lat + dlat
#             lon = base_lon + dlon
#             nodes.append({"id": i, "lat": lat, "lon": lon})
#         return nodes
        
#     nodes_data = generate_random_nodes(num_nodes, base_lat, base_lon, radius_km)
    
#     # 2. Queries Generation
#     events = []
#     for eid in range(num_events):
#         num_orders = random.randint(min_orders, max_orders)
#         num_drivers = random.randint(1, max_drivers)
#         depot_node = random.randrange(num_nodes)

#         orders = []
#         for oid in range(num_orders):
#             pickup = random.randrange(num_nodes)
#             dropoff = random.randrange(num_nodes)
#             while dropoff == pickup:
#                 dropoff = random.randrange(num_nodes)

#             price = round(random.uniform(100.0, 1000.0), 2)
#             is_premium = random.random() < 0.25
#             is_cancelled = random.random() < 0.05

#             orders.append({
#                 "order_id": oid,
#                 "pickup": pickup,
#                 "dropoff": dropoff,
#                 "order_price": price,
#                 "is_cancelled": is_cancelled,
#                 "is_premium": is_premium
#             })

#         events.append({
#             "event_id": eid,
#             "orders": orders,
#             "fleet": {
#                 "num_delievery_guys": num_drivers,
#                 "depot_node": depot_node
#             }
#         })

#     queries_json_content = {
#         "meta": {
#             "num_events": num_events
#         },
#         "events": events
#     }
    
#     # --- Process Queries using C++ logic simulation ---
    
#     city_graph = Graph(nodes_data)
    
#     results = []
#     for event_data in queries_json_content["events"]:
#         # Simulate C++ processing time
#         processing_time = round(random.uniform(50.0, 500.0), 3) 
        
#         # Simulate C++ scheduler.process_query
#         result = process_query_in_python(event_data, city_graph)
#         result["processing_time"] = processing_time # Add simulated time
#         results.append(result)

#     # --- Final Output Structure (Matching C++ main.cpp) ---
#     final_output = {
#         "meta": queries_json_content["meta"],
#         "results": results
#     }

#     # Use json.dumps for the final output string
#     return json.dumps(final_output, indent=4)

# print(generate_expected_output_json())










import json
import math
import random

# Configuration
N_NODES = 100
N_EVENTS = 3
MIN_ORDERS_PER_EVENT = 3
MAX_ORDERS_PER_EVENT = 7
DEPOT_NODE = 0

BASE_LAT = 19.0760   # roughly Mumbai, just to look "real"
BASE_LON = 72.8777

random.seed(42)

def generate_nodes():
    nodes = []
    for i in range(N_NODES):
        # random jitter of up to ~3km in each direction
        dlat = random.uniform(-0.03, 0.03)
        dlon = random.uniform(-0.03, 0.03)
        nodes.append({
            "id": i,
            "lat": BASE_LAT + dlat,
            "lon": BASE_LON + dlon
        })
    return nodes

def approx_length_m(a, b):
    # cheap distance approximation in meters using lat/lon
    dlat = (a["lat"] - b["lat"]) * 111_000.0
    # longitude length shrinks by cos(lat); use base lat
    dlon = (a["lon"] - b["lon"]) * 111_000.0 * math.cos(math.radians(BASE_LAT))
    return math.hypot(dlat, dlon)

def generate_edges(nodes, k_nearest=4):
    edges = []
    edge_id = 0
    N = len(nodes)
    for i in range(N):
        dists = []
        for j in range(N):
            if i == j:
                continue
            d = approx_length_m(nodes[i], nodes[j])
            dists.append((d, j))
        dists.sort()
        # connect to k nearest neighbors
        for d, j in dists[:k_nearest]:
            if i < j:  # avoid duplicates
                length = d
                # random average speed between 5 and 15 m/s
                speed = random.uniform(5.0, 15.0)
                avg_time = length / speed
                edges.append({
                    "id": edge_id,
                    "u": i,
                    "v": j,
                    "length": length,
                    "average_time": avg_time,
                    "oneway": False,
                    "type": "road"
                })
                edge_id += 1
    return edges

def generate_graph():
    nodes = generate_nodes()
    edges = generate_edges(nodes)
    graph = {
        "meta": {
            "nodes": len(nodes)
        },
        "nodes": nodes,
        "edges": edges
    }
    with open("graph.json", "w") as f:
        json.dump(graph, f, indent=2)
    print(f"graph.json written with {len(nodes)} nodes and {len(edges)} edges")

def generate_queries():
    events = []
    global_order_id = 1
    for _ in range(N_EVENTS):
        n_orders = random.randint(MIN_ORDERS_PER_EVENT, MAX_ORDERS_PER_EVENT)
        orders = []
        for _ in range(n_orders):
            pickup = random.randint(1, N_NODES - 1)
            dropoff = random.randint(1, N_NODES - 1)
            while dropoff == pickup:
                dropoff = random.randint(1, N_NODES - 1)
            price = random.uniform(150.0, 800.0)
            is_premium = random.random() < 0.3   # ~30% premium
            is_cancelled = random.random() < 0.1 # ~10% cancelled

            orders.append({
                "order_id": global_order_id,
                "pickup": pickup,
                "dropoff": dropoff,
                "order_price": round(price, 2),
                "is_premium": is_premium,
                "is_cancelled": is_cancelled
            })
            global_order_id += 1

        fleet = {
            "num_delievery_guys": random.randint(1, 4),
            "depot_node": DEPOT_NODE
        }

        events.append({
            "orders": orders,
            "fleet": fleet
        })

    queries = {
        "meta": {
            "description": "Randomly generated test events",
            "num_events": len(events)
        },
        "events": events
    }

    with open("queries.json", "w") as f:
        json.dump(queries, f, indent=2)
    print(f"queries.json written with {len(events)} events")

def generate_expected_output():
    with open("queries.json") as f:
        queries = json.load(f)

    results = []
    for event in queries["events"]:
        # ignore cancelled orders in this baseline
        orders = [o for o in event["orders"] if not o["is_cancelled"]]
        depot = event["fleet"]["depot_node"]
        route = [depot]
        order_ids = []

        # Very simple baseline: single driver, sequential pickups and dropoffs
        for o in orders:
            route.append(o["pickup"])
            route.append(o["dropoff"])
            order_ids.append(o["order_id"])

        # We don't know the real shortest path here; just use 0 as placeholder
        total_delivery_time = 0.0

        result = {
            "assignments": [
                {
                    "driver_id": 0,
                    "route": route,
                    "order_ids": order_ids
                }
            ],
            "metrics": {
                "total_delivery_time_s": total_delivery_time
            },
            "processing_time": 0.0
        }
        results.append(result)

    expected = {
        "results": results
    }

    with open("expected_output.json", "w") as f:
        json.dump(expected, f, indent=2)
    print(f"expected_output.json written with {len(results)} results")

if __name__ == "__main__":
    generate_graph()
    generate_queries()
    generate_expected_output()

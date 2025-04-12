const apiUrl = "http://localhost:8080/api/";

let startId = -1;
let endId = -1;
let pathGroup = undefined;
let map = undefined;
let polylineLayer = undefined;

function mapInit() {
    // coordinates of the Computer Science Building of Uni Stuttgart
    map = L.map("map").setView([48.7451545063299, 9.106633428929493], 17);

    pathGroup = L.layerGroup().addTo(map);

    L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
        attribution:
            '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
    }).addTo(map);

    var marker1 = L.marker([48.7451545063299, 9.106633428929493], {
        draggable: true,
    }).addTo(map);
    var marker2 = L.marker([48.7451545063299, 9.106633428929493], {
        draggable: true,
    }).addTo(map);

    function updateCoordinates() {
        // Get the current positions of the markers
        var latlng1 = marker1.getLatLng();
        var latlng2 = marker2.getLatLng();

        // get nearest node from the graph
        fetch(
            apiUrl + "nearest_node?lat=" + latlng1.lat + "&lon=" + latlng1.lng
        )
            .then((response) => {
                if (response.ok) {
                    response.json().then((data) => {
                        console.log(data);
                        startId = data.node;
                        // Update the marker position to the nearest node
                        marker1.setLatLng([data.lat, data.lon]);
                        // Update the HTML content with the nearest node coordinates
                        document.getElementById("start-coords").textContent =
                            data.lat + ", " + data.lon;
                    });
                } else {
                    console.error(
                        "Error fetching nearest node data:",
                        response.statusText
                    );
                }
            })
            .catch((error) => {
                console.error("Error:", error);
            });

        fetch(
            apiUrl + "nearest_node?lat=" + latlng2.lat + "&lon=" + latlng2.lng
        )
            .then((response) => {
                if (response.ok) {
                    response.json().then((data) => {
                        console.log(data);
                        endId = data.node;
                        // Update the marker position to the nearest node
                        marker2.setLatLng([data.lat, data.lon]);
                        // Update the HTML content with the nearest node coordinates
                        document.getElementById("end-coords").textContent =
                            data.lat + ", " + data.lon;
                    });
                } else {
                    console.error(
                        "Error fetching nearest node data:",
                        response.statusText
                    );
                }
            })
            .catch((error) => {
                console.error("Error:", error);
            });
    }

    // Initially update coordinates on page load
    updateCoordinates();

    // Add event listeners to update coordinates when a marker is dragged
    marker1.on("dragend", updateCoordinates);
    marker2.on("dragend", updateCoordinates);
}

function init() {
    fetch(apiUrl).then((response) => {
        if (response.ok) {
            response.json().then((data) => {
                const type = data.type;
                if (type == "simple") {
                    // only show the dijkstra and clear button
                    document.getElementById("complex-ch").style.display =
                        "none";
                    document.getElementById("complex-hub-label").style.display =
                        "none";
                } else {
                    // dont show the dijkstra button
                    document.getElementById("simple-dijkstra").style.display =
                        "none";
                }
            });
        } else {
            console.error("Error fetching data:", response.statusText);
        }
    });
}

function dijkstraQuery() {
    if (startId == -1 || endId == -1) {
        alert("Please select start and end coordinates!");
        return;
    }

    disableButtons();

    fetch(apiUrl + "dijkstra?start=" + startId + "&end=" + endId).then(
        (response) => {
            if (response.ok) {
                response.json().then((data) => {
                    console.log(data);
                    const path = data.path;
                    const distance = data.distance;
                    console.log("Distance: " + distance);

                    clearPath();
                    polylineLayer = L.polyline(path.coordinates).addTo(map);
                    // Adjust the map view to show the complete polyline
                    map.fitBounds(polylineLayer.getBounds());
                });
            } else {
                console.error(
                    "Error fetching dijkstra data:",
                    response.statusText
                );
            }
        }
    );

    enableButtons();
}

function chQuery() {
    if (startId == -1 || endId == -1) {
        alert("Please select start and end coordinates!");
        return;
    }

    disableButtons();

    fetch(apiUrl + "ch?start=" + startId + "&end=" + endId).then((response) => {
        if (response.ok) {
            response.json().then((data) => {
                console.log(data);
                const path = data.path;
                const distance = data.distance;
                console.log("Distance: " + distance);

                clearPath();
                polylineLayer = L.polyline(path.coordinates).addTo(map);
                // Adjust the map view to show the complete polyline
                map.fitBounds(polylineLayer.getBounds());
                enableButtons();
            });
        } else {
            enableButtons();
            console.error("Error fetching dijkstra data:", response.statusText);
        }
    });
}

function clearPath() {
    if (polylineLayer != undefined && map.hasLayer(polylineLayer))
        map.removeLayer(polylineLayer);
}

function graphRequest() {
    // Disable all buttons to prevent multiple submissions
    disableButtons();
    fetch(apiUrl).then((response) => {
        if (response.ok) {
            console.log(response);
            response.json().then((data) => {
                console.log(data);
                console.log(data["status"]);
            });
        } else {
            console.error("Error fetching graph data:", response.statusText);
        }
    });

    // Re-enable buttons after the request is complete
    enableButtons();
}

// Helper function to disable all buttons on the page
function disableButtons() {
    document.querySelectorAll("button").forEach((btn) => {
        btn.disabled = true;
        // Optionally add a visual indicator using CSS classes or inline styles
        btn.classList.add("opacity-50", "cursor-not-allowed");
    });
}

// Helper function to re-enable all buttons on the page
function enableButtons() {
    document.querySelectorAll("button").forEach((btn) => {
        btn.disabled = false;
        btn.classList.remove("opacity-50", "cursor-not-allowed");
    });
}

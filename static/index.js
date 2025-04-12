const apiUrl = "http://localhost:8080/api/graph";

function mapInit() {
    // coordinates of the Computer Science Building of Uni Stuttgart
    var map = L.map("map").setView([48.7451545063299, 9.106633428929493], 17);

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

        // Format the coordinates to 5 decimal places and update the HTML content
        document.getElementById("marker1-coords").textContent =
            latlng1.lat.toFixed(5) + ", " + latlng1.lng.toFixed(5);
        document.getElementById("marker2-coords").textContent =
            latlng2.lat.toFixed(5) + ", " + latlng2.lng.toFixed(5);
    }

    // Initially update coordinates on page load
    updateCoordinates();

    // Add event listeners to update coordinates when a marker is dragged
    marker1.on("dragend", updateCoordinates);
    marker2.on("dragend", updateCoordinates);
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

    setTimeout(enableButtons, 5000);
    // Re-enable buttons after the request is complete
    //enableButtons();
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

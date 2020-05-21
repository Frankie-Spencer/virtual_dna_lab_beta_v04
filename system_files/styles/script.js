

$(".card_agent").hover(function (event) {

    var currentValue = $(this).children(".agent_text").html();

    // var parentElement = $(this).closest(".card_complex_group");

    $(this).closest(".card_complex_group").find(".card_agent").filter(function () {

        var child = $(this).children(".agent_text");

        if (child.length > 0) {

            return child[0].innerHTML === currentValue

        }

        return false;

    }).css("background-color", event.type === "mouseenter" ? "green" : "white");

});


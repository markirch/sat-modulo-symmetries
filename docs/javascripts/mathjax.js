window.MathJax = {
    tex: {
        inlineMath: [["\\(", "\\)"]],
        displayMath: [["\\[", "\\]"]],
        processEscapes: true,
        processEnvironments: true
    },
    options: {
        ignoreHtmlClass: ".*|",
        processHtmlClass: "arithmatex"
    },
    chtml: {
        scale: 0.85
    },
    svg: {
        scale: 0.85
    }
};

document$.subscribe(() => {
    MathJax.typesetPromise()
})
// KaTeX-ініціалізація для MkDocs-Material (pymdownx.arithmatex: generic).
// Рендер після кожної навігації (Material вантажить сторінки через instant-loading).
document$.subscribe(function () {
  renderMathInElement(document.body, {
    delimiters: [
      { left: "$$",  right: "$$",  display: true  },
      { left: "$",   right: "$",   display: false },
      { left: "\\(", right: "\\)", display: false },
      { left: "\\[", right: "\\]", display: true  }
    ],
    // формули тому щільні, але прості за класом: не валимо сторінку на промаху
    throwOnError: false
  });
});

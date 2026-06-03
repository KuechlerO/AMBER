(function () {
    var STORAGE_KEY = 'amber-theme';
    var DEFAULT_THEME = 'amber';
    var VALID = ['green', 'amber', 'scientific'];

    function applyTheme(theme) {
        if (VALID.indexOf(theme) === -1) {
            theme = DEFAULT_THEME;
        }
        document.documentElement.setAttribute('data-theme', theme);
        localStorage.setItem(STORAGE_KEY, theme);
        document.querySelectorAll('.theme-option').forEach(function (btn) {
            var active = btn.getAttribute('data-theme') === theme;
            btn.classList.toggle('theme-option-active', active);
            btn.setAttribute('aria-pressed', active ? 'true' : 'false');
        });
    }

    document.addEventListener('DOMContentLoaded', function () {
        document.querySelectorAll('.theme-option').forEach(function (btn) {
            btn.addEventListener('click', function () {
                applyTheme(btn.getAttribute('data-theme'));
            });
        });
        var saved = localStorage.getItem(STORAGE_KEY) || DEFAULT_THEME;
        applyTheme(saved);
    });

    window.amberApplyTheme = applyTheme;
})();

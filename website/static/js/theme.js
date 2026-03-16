// Theme Switcher - Light/Dark Mode with Smooth Transitions

// Get saved theme or default to light
let currentTheme = localStorage.getItem('theme') || 'light';

// Apply theme on page load
document.addEventListener('DOMContentLoaded', function() {
    applyTheme(currentTheme);
});

// Toggle theme function
function toggleTheme(theme) {
    currentTheme = theme;
    applyTheme(theme);
    localStorage.setItem('theme', theme);
}

// Apply theme to body and update buttons with smooth transition
function applyTheme(theme) {
    const body = document.body;
    const lightBtn = document.getElementById('theme-light');
    const darkBtn = document.getElementById('theme-dark');
    
    // Add transition class before changing theme
    body.style.transition = 'all 0.5s cubic-bezier(0.4, 0, 0.2, 1)';
    
    if (theme === 'dark') {
        body.classList.add('dark-mode');
        if (lightBtn) lightBtn.classList.remove('active');
        if (darkBtn) darkBtn.classList.add('active');
    } else {
        body.classList.remove('dark-mode');
        if (darkBtn) darkBtn.classList.remove('active');
        if (lightBtn) lightBtn.classList.add('active');
    }
    
    // Trigger custom event for chart updates with slight delay for smooth transition
    setTimeout(() => {
        window.dispatchEvent(new CustomEvent('themeChanged', { detail: { theme: theme } }));
    }, 100);
}

// Get current theme colors for charts
function getChartColors() {
    const isDark = document.body.classList.contains('dark-mode');
    
    return {
        background: isDark ? '#1e293b' : '#ffffff',
        gridColor: isDark ? '#334155' : '#e2e8f0',
        textColor: isDark ? '#f1f5f9' : '#0f172a',
        primaryColor: isDark ? '#3b82f6' : '#2563eb'
    };
}

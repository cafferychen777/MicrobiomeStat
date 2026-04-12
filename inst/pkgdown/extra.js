// Enhanced pkgdown functionality for MicrobiomeStat
document.addEventListener('DOMContentLoaded', function() {
  
  // Smooth scrolling for anchor links
  document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
      e.preventDefault();
      const target = document.querySelector(this.getAttribute('href'));
      if (target) {
        target.scrollIntoView({
          behavior: 'smooth',
          block: 'start'
        });
      }
    });
  });
  
  // Add navbar scroll effect
  let lastScroll = 0;
  const navbar = document.querySelector('.navbar');
  
  if (navbar) {
    navbar.style.transition = 'transform 0.3s ease, box-shadow 0.3s ease';
    
    window.addEventListener('scroll', () => {
      const currentScroll = window.pageYOffset;
      
      if (currentScroll <= 0) {
        navbar.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
        return;
      }
      
      if (currentScroll > lastScroll && currentScroll > 100) {
        // Scrolling down
        navbar.style.transform = 'translateY(-100%)';
      } else {
        // Scrolling up
        navbar.style.transform = 'translateY(0)';
        navbar.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
      }
      
      lastScroll = currentScroll;
    });
  }
  
  // Add copy button to code blocks
  const codeBlocks = document.querySelectorAll('pre');
  
  codeBlocks.forEach(block => {
    // Skip if button already exists
    if (block.querySelector('.copy-btn')) return;
    
    const button = document.createElement('button');
    button.className = 'btn btn-sm btn-outline-secondary copy-btn';
    button.textContent = 'Copy';
    button.style.cssText = 'position: absolute; top: 0.5rem; right: 0.5rem; opacity: 0; transition: opacity 0.3s; z-index: 10;';
    
    block.style.position = 'relative';
    block.appendChild(button);
    
    block.addEventListener('mouseenter', () => {
      button.style.opacity = '1';
    });
    
    block.addEventListener('mouseleave', () => {
      button.style.opacity = '0';
    });
    
    button.addEventListener('click', () => {
      const code = block.querySelector('code');
      if (!code) return;
      
      const text = code.textContent || code.innerText;
      
      if (navigator.clipboard) {
        navigator.clipboard.writeText(text).then(() => {
          button.textContent = 'Copied!';
          button.classList.remove('btn-outline-secondary');
          button.classList.add('btn-success');
          
          setTimeout(() => {
            button.textContent = 'Copy';
            button.classList.remove('btn-success');
            button.classList.add('btn-outline-secondary');
          }, 2000);
        }).catch(err => {
          console.error('Failed to copy:', err);
        });
      } else {
        // Fallback for older browsers
        const textarea = document.createElement('textarea');
        textarea.value = text;
        textarea.style.position = 'fixed';
        textarea.style.left = '-999999px';
        document.body.appendChild(textarea);
        textarea.select();
        
        try {
          document.execCommand('copy');
          button.textContent = 'Copied!';
          button.classList.remove('btn-outline-secondary');
          button.classList.add('btn-success');
          
          setTimeout(() => {
            button.textContent = 'Copy';
            button.classList.remove('btn-success');
            button.classList.add('btn-outline-secondary');
          }, 2000);
        } catch (err) {
          console.error('Failed to copy:', err);
        }
        
        document.body.removeChild(textarea);
      }
    });
  });
  
  // Add fade-in animation to page elements
  const animateElements = document.querySelectorAll('.card, .ref-description, h2, h3');
  
  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        entry.target.classList.add('fade-in-up');
        observer.unobserve(entry.target);
      }
    });
  }, {
    threshold: 0.1,
    rootMargin: '0px 0px -50px 0px'
  });
  
  animateElements.forEach(element => {
    observer.observe(element);
  });
  
  // Enhance tables with responsive wrapper
  const tables = document.querySelectorAll('table');
  tables.forEach(table => {
    if (!table.parentElement.classList.contains('table-responsive')) {
      const wrapper = document.createElement('div');
      wrapper.className = 'table-responsive';
      table.parentNode.insertBefore(wrapper, table);
      wrapper.appendChild(table);
    }
  });
  
  // Add active state to TOC
  const tocLinks = document.querySelectorAll('#toc a');
  const sections = document.querySelectorAll('h2[id], h3[id]');
  
  if (tocLinks.length > 0 && sections.length > 0) {
    window.addEventListener('scroll', () => {
      let current = '';
      
      sections.forEach(section => {
        const sectionTop = section.offsetTop;
        const sectionHeight = section.clientHeight;
        if (window.pageYOffset >= sectionTop - 200) {
          current = section.getAttribute('id');
        }
      });
      
      tocLinks.forEach(link => {
        link.classList.remove('active');
        if (link.getAttribute('href') === '#' + current) {
          link.classList.add('active');
        }
      });
    });
  }
  
  // Improve search box
  const searchBox = document.querySelector('input[type="search"]');
  if (searchBox) {
    searchBox.setAttribute('placeholder', '🔍 Search functions, datasets, or topics...');
  }
  
  // Add keyboard shortcuts
  document.addEventListener('keydown', function(e) {
    // Ctrl/Cmd + K for search
    if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
      e.preventDefault();
      const searchBox = document.querySelector('input[type="search"]');
      if (searchBox) {
        searchBox.focus();
        searchBox.select();
      }
    }
    
    // Escape to close search
    if (e.key === 'Escape') {
      const searchBox = document.querySelector('input[type="search"]');
      if (searchBox && document.activeElement === searchBox) {
        searchBox.blur();
      }
    }
  });
  
  // Add loading indicator for slow operations
  const links = document.querySelectorAll('a[href$=".html"]');
  links.forEach(link => {
    link.addEventListener('click', function() {
      if (this.href && !this.href.includes('#')) {
        document.body.style.cursor = 'progress';
      }
    });
  });
  
  // Console welcome message
  console.log('%c🦠 Welcome to MicrobiomeStat!', 
    'color: #4CAF50; font-size: 20px; font-weight: bold;');
  console.log('%cComprehensive Statistical and Visualization Methods for Microbiome Data', 
    'color: #558B2F; font-size: 14px;');
  console.log('📚 Documentation: https://www.microbiomestat.wiki');
  console.log('💬 Discord: https://discord.gg/U6vNamyk');
  console.log('🐛 Issues: https://github.com/cafferychen777/MicrobiomeStat/issues');
});

// Add version check notification (optional)
window.addEventListener('load', function() {
  // Check if there's a newer version available
  const currentVersion = document.querySelector('.version')?.textContent;
  if (currentVersion && window.location.hostname !== 'localhost') {
    fetch('https://api.github.com/repos/cafferychen777/MicrobiomeStat/releases/latest')
      .then(response => response.json())
      .then(data => {
        const latestVersion = data.tag_name;
        if (latestVersion && latestVersion !== currentVersion) {
          console.log(`📦 New version ${latestVersion} is available!`);
        }
      })
      .catch(err => {
        // Silently fail if GitHub API is unavailable
      });
  }
});
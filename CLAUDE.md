# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is an **al-folio** academic website - a Jekyll-based static site theme designed for researchers and academics. The site is deployed to GitHub Pages and includes publications, CV, blog posts, projects, and other academic content.

**Owner**: Seok-Jin Kang
**URL**: https://seok-jin1.github.io
**Theme**: al-folio (https://github.com/alshedivat/al-folio)

## Development Commands

### Local Development (Docker - Recommended)

```bash
# Pull and run the site locally
docker compose pull
docker compose up

# Site will be available at http://localhost:8080

# Use slim image (faster, smaller)
docker compose -f docker-compose-slim.yml up

# Rebuild docker image from scratch
docker compose up --build

# Debug docker issues
docker compose up -d
docker compose logs
docker compose exec -it jekyll /bin/bash
./bin/entry_point.sh
```

### Building the Site

```bash
# Set production environment and build
export JEKYLL_ENV=production
bundle exec jekyll build

# Build output is in _site/

# Purge unused CSS after build
npm install -g purgecss
purgecss -c purgecss.config.js
```

### Git Flow

```bash
# 변경된 파일 스테이징
git add <파일경로>         # 특정 파일
git add .                  # 전체 변경사항

# 커밋
git commit -m "커밋 메시지"

# 푸시 (GitHub Actions가 자동으로 배포)
git push origin main
```

### Code Quality and Formatting

```bash
# Format code with Prettier (includes Liquid templates)
npx prettier --write .

# The repo uses pre-commit hooks for formatting checks
```

## Architecture

### Jekyll Collections

The site uses Jekyll collections for organizing academic content:

- **`_bibliography/`** - BibTeX files for publications (primary: `papers.bib`)
- **`_posts/`** - Blog posts in markdown format
- **`_projects/`** - Project pages
- **`_news/`** - News items displayed on homepage
- **`_books/`** - Book reviews and reading lists
- **`_pages/`** - Static pages (about, publications, CV, etc.)

### Key Directories

- **`_layouts/`** - Liquid layout templates (`about`, `bib`, `cv`, `post`, `distill`, etc.)
- **`_includes/`** - Reusable Liquid partials (header, footer, figure, citation, etc.)
- **`_sass/`** - SCSS stylesheets
- **`_plugins/`** - Custom Jekyll plugins (Ruby)
- **`assets/`** - Static assets (images, PDFs, JavaScript, CSS)
- **`_site/`** - Generated static site (git-ignored, created on build)

### Custom Jekyll Plugins

Located in `_plugins/`, written in Ruby:

- **`cache-bust.rb`** - Asset cache busting
- **`external-posts.rb`** - Fetch external blog posts
- **`google-scholar-citations.rb`** - Fetch Google Scholar citation counts
- **`inspirehep-citations.rb`** - Fetch InspireHEP citation counts
- **`details.rb`** - Details/summary HTML elements
- **`file-exists.rb`** - Check file existence
- **`hide-custom-bibtex.rb`** - Filter BibTeX fields
- **`remove-accents.rb`** - String normalization

### Bibliography System

Uses **Jekyll-Scholar** plugin to generate publication pages from BibTeX:

- Primary bib file: `_bibliography/papers.bib`
- Configuration in `_config.yml` under `scholar:` section
- Author name: Kang, Seok-Jin
- Bibliography template: `_layouts/bib.liquid`
- Supports features like: abstracts, PDFs, DOIs, code links, slides, posters
- Publication badges: Altmetric, Dimensions, Google Scholar, InspireHEP
- Custom fields: `pub_category` (featured/published), `preview` (thumbnail), `selected`, `impact_factor`, `featured_rank`

### Build Process

The GitHub Actions workflow (`.github/workflows/deploy.yml`) performs:

1. Checkout code
2. Setup Ruby 3.3.5 with bundler cache
3. Setup Python 3.13 for nbconvert (Jupyter notebooks)
4. Install ImageMagick (for responsive image processing)
5. Install Python dependencies (`nbconvert`)
6. Build Jekyll site with `JEKYLL_ENV=production`
7. Run PurgeCSS to remove unused CSS
8. Deploy to `gh-pages` branch

### Dependencies

**Ruby (Gemfile):**

- Jekyll with numerous plugins (scholar, feed, imagemagick, jupyter-notebook, etc.)
- Key plugins: `jekyll-scholar`, `jekyll-imagemagick`, `jekyll-jupyter-notebook`

**Node.js (package.json):**

- Prettier with Liquid plugin for code formatting
- Very minimal Node dependencies

**System:**

- ImageMagick (for image processing)
- Python with nbconvert (for Jupyter notebook support)

## Working with Publications

### Adding a New Publication

1. Edit `_bibliography/papers.bib`
2. Add BibTeX entry with these recommended fields:
   - Standard: `author`, `title`, `journal`, `year`, `volume`, `pages`, `doi`, `url`
   - Custom: `abbr` (journal abbreviation), `preview` (image filename), `selected` (true/false)
   - Optional: `pdf`, `code`, `slides`, `poster`, `abstract`, `arxiv`, `bibtex_show`
   - Categories: `pub_category` (featured/published), `featured_rank` (number)
3. Add preview image to `assets/img/publication_preview/` if using `preview` field
4. Add PDF to `assets/pdf/` if using `pdf` field
5. Jekyll-Scholar automatically regenerates publication pages on build

### Publication Display Logic

- Publications sorted by year (descending) by default
- `selected: true` makes publication appear on homepage
- `pub_category: featured` with `featured_rank` controls featured publications order
- Scholar config: `_config.yml` lines 266-330

## Site Configuration

Primary config file: **`_config.yml`**

Key settings to modify:

- Lines 5-8: Personal info (title, name)
- Line 21: Site URL
- Lines 266-268: Scholar author name (for publication highlighting)
- Lines 293-297: Publication badge toggles
- Lines 381-395: Feature flags (analytics, math, dark mode, etc.)

## Important Notes

- **Do not modify `_site/` directory** - it's auto-generated
- **Liquid templates** use `.liquid` extension (not `.html`)
- **Image optimization** is automatic via jekyll-imagemagick (creates responsive WebP images)
- **Math support** via MathJax (enable with `enable_math: true` in config)
- **Dark mode** is built-in and respects user OS preferences
- **GitHub Actions** auto-deploys on push to main branch
- The site requires **ImageMagick** installed for image processing to work

## File Naming Conventions

- Blog posts: `_posts/YYYY-MM-DD-title.md`
- Projects: `_projects/project-name.md`
- News: `_news/announcement-N.md`
- Pages: `_pages/pagename.md` with `permalink: /pagename/` in frontmatter

## Responsive Images

The jekyll-imagemagick plugin automatically generates responsive images:

- Input formats: JPG, JPEG, PNG, TIFF, GIF
- Output: WebP at multiple widths (480px, 800px, 1400px)
- Processes files in `assets/img/`
- Lazy loading enabled by default (`lazy_loading_images: true`)

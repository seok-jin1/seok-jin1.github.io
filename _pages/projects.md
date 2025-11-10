---
layout: page
title: Projects
permalink: /projects/
description: Research projects and funded work.
nav: true
nav_order: 5
---

{% assign projects = site.data.research_projects %}

<div class="research-projects">
  {% for project in projects %}
    <div class="project-item">
      <div class="project-header">
        <h3 class="project-title">{{ project.title }}</h3>
        <div class="project-meta">
          <span class="project-institution">{{ project.institution }}</span>
          <span class="project-year">{{ project.year }}</span>
        </div>
      </div>
      {% if project.description %}
        <ul class="project-description">
          {% for item in project.description %}
            <li>{{ item }}</li>
          {% endfor %}
        </ul>
      {% endif %}
    </div>
  {% endfor %}
</div>

<style>
.research-projects {
  margin-top: 2rem;
}

.project-item {
  margin-bottom: 2.5rem;
  padding-bottom: 1.5rem;
  border-bottom: 1px solid var(--global-divider-color);
}

.project-item:last-child {
  border-bottom: none;
}

.project-header {
  margin-bottom: 0.75rem;
}

.project-title {
  font-size: 1.25rem;
  font-weight: 600;
  margin-bottom: 0.5rem;
  color: var(--global-text-color);
}

.project-meta {
  display: flex;
  gap: 1rem;
  flex-wrap: wrap;
  font-size: 0.95rem;
  color: var(--global-text-color-light);
}

.project-institution {
  font-weight: 500;
}

.project-year {
  color: var(--global-theme-color);
}

.project-description {
  margin-top: 0.5rem;
  margin-left: 1.25rem;
  line-height: 1.6;
}

.project-description li {
  margin-bottom: 0.25rem;
}
</style>

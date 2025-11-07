---
layout: page
permalink: /books/
title: Books
description: curated book publications organized by international and domestic releases.
nav: true
nav_order: 3
---

{% assign books = site.data.books %}

<div class="books">
  {% if books.international %}
    <section>
      <h2>International Publications</h2>
      <ul class="book-list">
        {% for book in books.international %}
          <li class="book-card z-depth-1">
            {% if book.cover %}
              <div class="book-cover">
                <img src="{{ book.cover | relative_url }}" alt="Cover of {{ book.title }}">
              </div>
            {% endif %}
            <div class="book-details">
              <div class="book-header">
                <div>
                  <div class="book-title">{{ book.title }}</div>
                  {% if book.subtitle %}
                    <div class="book-subtitle">{{ book.subtitle }}</div>
                  {% endif %}
                </div>
                {% if book.publication_date %}
                  <div class="book-date">{{ book.publication_date }}</div>
                {% endif %}
              </div>
              <div class="book-body">
                {% if book.publisher %}
                  <div class="book-meta"><strong>Publisher:</strong> {{ book.publisher }}</div>
                {% endif %}
                {% if book.isbn %}
                  <div class="book-meta"><strong>ISBN:</strong> {{ book.isbn }}</div>
                {% endif %}
              </div>
              {% if book.chapters %}
                <div class="book-chapters">
                  <div class="book-meta"><strong>Chapter authors</strong></div>
                  <ul>
                    {% for chapter in book.chapters %}
                      {% assign highlight = false %}
                      {% if chapter.author contains 'Seok-Jin Kang' or chapter.author contains '강석진' %}
                        {% assign highlight = true %}
                      {% endif %}
                      <li>
                        <span class="chapter-author">
                          {% if highlight %}<strong>{{ chapter.author }}</strong>{% else %}{{ chapter.author }}{% endif %}
                        </span>
                        <span class="chapter-topic">
                          — {% if highlight %}<strong>{{ chapter.topic }}</strong>{% else %}{{ chapter.topic }}{% endif %}
                        </span>
                        {% if chapter.author_en or chapter.topic_en %}
                          <span class="chapter-translation">
                            (
                            {% if chapter.author_en %}
                              {% if highlight %}<strong>{{ chapter.author_en }}</strong>{% else %}{{ chapter.author_en }}{% endif %}
                            {% endif %}
                            {% if chapter.topic_en %}
                              {% if chapter.author_en %} — {% endif %}
                              {% if highlight %}<strong>{{ chapter.topic_en }}</strong>{% else %}{{ chapter.topic_en }}{% endif %}
                            {% endif %}
                            )
                          </span>
                        {% endif %}
                      </li>
                    {% endfor %}
                  </ul>
                </div>
              {% endif %}
              <div class="book-actions">
                {% if book.url %}
                  <a class="btn btn-sm z-depth-0" href="{{ book.url }}" target="_blank" rel="noopener">View book</a>
                {% else %}
                  <span class="book-url-pending">URL coming soon</span>
                {% endif %}
              </div>
            </div>
          </li>
        {% endfor %}
      </ul>
    </section>
  {% endif %}

{% if books.domestic %}

<section>
<h2>Domestic Publications</h2>
<ul class="book-list">
{% for book in books.domestic %}
<li class="book-card z-depth-1">
{% if book.cover %}
<div class="book-cover">
<img src="{{ book.cover | relative_url }}" alt="{{ book.title }} 표지">
</div>
{% endif %}
<div class="book-details">
<div class="book-header">
<div>
<div class="book-title">{{ book.title }}</div>
{% if book.subtitle %}
<div class="book-subtitle">{{ book.subtitle }}</div>
{% endif %}
</div>
{% if book.publication_date %}
<div class="book-date">{{ book.publication_date }}</div>
{% endif %}
</div>
<div class="book-body">
{% if book.publisher %}
<div class="book-meta"><strong>Publisher:</strong> {{ book.publisher }}</div>
{% endif %}
{% if book.isbn %}
<div class="book-meta"><strong>ISBN:</strong> {{ book.isbn }}</div>
{% endif %}
</div>
{% if book.chapters %}
<div class="book-chapters">
<div class="book-meta"><strong>Chapter authors</strong></div>
<ul>
{% for chapter in book.chapters %}
{% assign highlight = false %}
{% if chapter.author contains 'Seok-Jin Kang' or chapter.author contains '강석진' %}
{% assign highlight = true %}
{% endif %}
<li>
<span class="chapter-author">
{% if highlight %}<strong>{{ chapter.author }}</strong>{% else %}{{ chapter.author }}{% endif %}
</span>
<span class="chapter-topic">
— {% if highlight %}<strong>{{ chapter.topic }}</strong>{% else %}{{ chapter.topic }}{% endif %}
</span>
{% if chapter.author_en or chapter.topic_en %}
<span class="chapter-translation">
(
{% if chapter.author_en %}
{% if highlight %}<strong>{{ chapter.author_en }}</strong>{% else %}{{ chapter.author_en }}{% endif %}
{% endif %}
{% if chapter.topic_en %}
{% if chapter.author_en %} — {% endif %}
{% if highlight %}<strong>{{ chapter.topic_en }}</strong>{% else %}{{ chapter.topic_en }}{% endif %}
{% endif %}
)
</span>
{% endif %}
</li>
{% endfor %}
</ul>
</div>
{% endif %}
<div class="book-actions">
{% if book.url %}
<a class="btn btn-sm z-depth-0" href="{{ book.url }}" target="_blank" rel="noopener">View book</a>
{% else %}
<span class="book-url-pending">URL coming soon</span>
{% endif %}
</div>
</div>
</li>
{% endfor %}
</ul>
</section>
{% endif %}

</div>

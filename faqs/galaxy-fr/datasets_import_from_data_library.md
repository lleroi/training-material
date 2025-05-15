---
title: Importing data from a data library
area: data upload
box_type: tip
layout: faq
contributors: [bebatut,shiltemann,nsoranzo,hexylena,wm75,paulineauffret]
---


As an alternative to uploading the data from a URL or your computer, the files may also have been made available from a *shared data library*:

1. Go into **Data (_Données partagées_)** (top panel) then **Data libraries (_Bibliothèque de données_)**
2. Navigate to **{{ include.path | default: "the correct folder as indicated by your instructor." }}** 
   {% unless include.path %}- On most Galaxies tutorial data will be provided in a folder named **GTN - Material --> Topic Name -> Tutorial Name**. {% endunless %}
3. Select **{{ include.selection | default: "the desired files" }}**
4. Click on **Export to History** {% icon galaxy-dropdown %} near the top and select **{{ include.astype | default: "as Datasets" }}** from the dropdown menu
5. In the pop-up window, choose
   {% if include.collection_type %}
   * *"Collection type"*: **{{ include.collection_type }}**
   {% endif %}
   * *"Select history"*: {% if include.tohistory %}{{ include.tohistory }}{% else %}the history you want to import the data to (or create a new one){% endif %}
6. Click on {% if include.collection_type %}**Continue**{% else %}**Import**{% endif %}
{% if include.collection_type %}
7. In the next dialog, give a suitable **Name**{% if include.collection_name %}, like `{{ include.collection_name }}`,{% endif %} to the new collection
8. Click on **Create collection**
{% endif %}

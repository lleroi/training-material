---
title: Creating a dataset collection
area: collections
box_type: tip
layout: faq
contributors: [shiltemann, hexylena, paulineauffret]
---

* Click on {% icon galaxy-selector %} **Select Items (_Opérer sur plusieurs jeux de données en même temps_)** at the top of the history panel
![Select Items button]({% link topics/galaxy-interface/images/historyItemControls.png %})
* Select {% if include.datasets_description %}_**{{ include.datasets_description }}**_{% else %}Check all the datasets in your history you would like to include{% endif %}
* Click on **Pour toute la sélection** {% icon galaxy-dropdown %}
* Choose **Build Dataset List**

  ![build list collection menu item]({% link topics/galaxy-interface/images/buildList.png %}){:width="15%"}

* Enter the name for your collection {% if include.name %} : **{{ include.name }}** {% endif %}
* Click **Create list** to build your collection
* Click on the checkmark icon at the top of your history again


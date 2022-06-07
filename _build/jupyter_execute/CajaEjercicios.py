#!/usr/bin/env python
# coding: utf-8

# # Ejercicios partícula en una caja

# In[1]:


from OptMultiple import MultipleChoice


# In[2]:


question = "La energía del estado base de una partícula en una caja es cero."
answers = [
    "Falso",
    "Cierto"
]
explanation = (
    "Dado que si E=0 entonces ψ=0"
    " y dicha solución no cumple que su cuadrado sea una densidad de probabilidad."
)
MultipleChoice(
    question, answers, correct_answer=0, explanation=explanation
)


# In[ ]:





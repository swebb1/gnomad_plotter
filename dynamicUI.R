
## VdVp


## Dynamically add annotations
annoCount <- reactiveVal(0)
annoRM <- reactiveVal(vector())

observeEvent(input$addAnno,{
  
  #btn=sum(input$addAnno, 1)
  annoCount(annoCount()+1)
  btn=annoCount()
  
  insertUI(
    selector = '#placeholder',
    ui = div(
      id = paste0("div_",btn),
      h5(paste0("Annotation ",btn)),
      textInput(paste0("annoName_",btn),label = "Name"),
      textAreaInput(paste0("annoPos_",btn),label = "Positions"),
      helpText("Add new positions on each line"),
      actionButton(paste0("remove_",btn),label = "Remove")
    )
  )
  
  observeEvent(input[[paste0("remove_",btn)]], {
    rmv = paste0("#div_",btn)
    removeUI(
      selector = rmv
    )
    annoRM(append(annoRM(),btn))
  })
  
})

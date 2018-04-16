export default {
  list: [
    {
      category: 'Sequence design',
      scenarios: [
        require('./SculptASequence.vue').default,
        require('./DesignOverhangs.vue').default,
        require('./SimulateGGAssemblies').default,
        require('./SwapDonorVectorPart').default,
        require('./SketchConstructs').default,
        require('./ConvertSequenceFiles').default
      ]
    },
    {
      category: 'Quality Control',
      scenarios: [
        require('./PredictDigestions').default,
        require('./SelectDigestions').default,
        require('./AnalyzeDigests').default,
        require('./SelectPrimers').default
      ]
    },
    {
      category: 'Sequence analysis',
      scenarios: [
        require('./PlotSequenceFeatures.vue').default,
        require('./EvaluateManufacturability.vue').default,
        require('./FindCommonBlocks.vue').default,
        require('./RenderSequenticons.vue').default
      ]
    }
  ]
}

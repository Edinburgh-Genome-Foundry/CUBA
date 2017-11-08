export default {
  list: [
    {
      category: 'Sequence design',
      scenarios: [
        require('./SculptASequence.vue'),
        require('./DesignOverhangs.vue'),
        require('./SimulateGGAssemblies'),
        require('./SwapDonorVectorPart'),
        require('./SketchConstructs')
      ]
    },
    {
      category: 'Quality Control',
      scenarios: [
        require('./PredictDigestions'),
        require('./SelectDigestions'),
        require('./AnalyzeDigests')
      ]
    },
    {
      category: 'Sequence analysis',
      scenarios: [
        require('./EvaluateManufacturability.vue'),
        require('./FindCommonBlocks.vue')
      ]
    }
  ]
}

export default {
  list: [
    {
      category: 'Sequence design',
      scenarios: [
        require('./SculptASequence.vue'),
        require('./DesignOverhangs.vue'),
        require('./SimulateGGAssemblies')
      ]
    },
    {
      category: 'Sequence analysis',
      scenarios: [
        require('./EvaluateManufacturability.vue'),
        require('./FindCommonBlocks.vue')
      ]
    },
    {
      category: 'Quality Control',
      scenarios: [
        require('./DigestionPatternPredictor'),
        require('./SelectDigestions')
      ]
    }
  ]
}
